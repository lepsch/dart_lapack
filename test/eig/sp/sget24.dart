      void sget24(COMP, JTYPE, THRESH, final Array<int> ISEED, NOUNIT, N, final Matrix<double> A, final int LDA, H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, final Matrix<double> VS, final int LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, RESULT, final Array<double> WORK, final int LWORK, IWORK, BWORK, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               COMP;
      int                INFO, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT;
      double               RCDEIN, RCDVIN, THRESH;
      bool               BWORK( * );
      int                ISEED( 4 ), ISLCT( * ), IWORK( * );
      double               A( LDA, * ), H( LDA, * ), HT( LDA, * ), RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ), WI( * ), WIT( * ), WITMP( * ), WORK( * ), WR( * ), WRT( * ), WRTMP( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               EPSIN;
      const              EPSIN = 5.9605e-8 ;
      String             SORT;
      int                I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, LIWORK, RSUB, SDIM, SDIM1;
      double               ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, SMLNUM, TMP, TOL, TOLIN, ULP, ULPINV, V, VIMIN, VRMIN, WNORM;
      int                IPNT( 20 );
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double               SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. External Functions ..
      //- bool               SSLECT;
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SSLECT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEESX, SGEMM, SLACPY, SORT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SIGN, SQRT

      // Check for errors

      INFO = 0;
      if ( THRESH < ZERO ) {
         INFO = -3;
      } else if ( NOUNIT <= 0 ) {
         INFO = -5;
      } else if ( N < 0 ) {
         INFO = -6;
      } else if ( LDA < 1 || LDA < N ) {
         INFO = -8;
      } else if ( LDVS < 1 || LDVS < N ) {
         INFO = -18;
      } else if ( LWORK < 3*N ) {
         INFO = -26;
      }

      if ( INFO != 0 ) {
         xerbla('SGET24', -INFO );
         return;
      }

      // Quick return if nothing to do

      for (I = 1; I <= 17; I++) { // 10
         RESULT[I] = -ONE;
      } // 10

      if (N == 0) return;

      // Important constants

      SMLNUM = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Precision' );
      ULPINV = ONE / ULP;

      // Perform tests (1)-(13)

      SELOPT = 0;
      LIWORK = N*N;
      for (ISORT = 0; ISORT <= 1; ISORT++) { // 120
         if ( ISORT == 0 ) {
            SORT = 'N';
            RSUB = 0;
         } else {
            SORT = 'S';
            RSUB = 6;
         }

         // Compute Schur form and Schur vectors, and test them

         slacpy('F', N, N, A, LDA, H, LDA );
         sgeesx('V', SORT, SSLECT, 'N', N, H, LDA, SDIM, WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[1+RSUB] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX1', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX1', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            return;
         }
         if ( ISORT == 0 ) {
            scopy(N, WR, 1, WRTMP, 1 );
            scopy(N, WI, 1, WITMP, 1 );
         }

         // Do Test (1) or Test (7)

         RESULT[1+RSUB] = ZERO;
         for (J = 1; J <= N - 2; J++) { // 30
            for (I = J + 2; I <= N; I++) { // 20
               if( H( I, J ) != ZERO ) RESULT( 1+RSUB ) = ULPINV;
            } // 20
         } // 30
         for (I = 1; I <= N - 2; I++) { // 40
            if( H( I+1, I ) != ZERO && H( I+2, I+1 ) != ZERO ) RESULT( 1+RSUB ) = ULPINV;
         } // 40
         for (I = 1; I <= N - 1; I++) { // 50
            if ( H( I+1, I ) != ZERO ) {
               if( H( I, I ) != H( I+1, I+1 ) || H( I, I+1 ) == ZERO || sign( ONE, H( I+1, I ) ) == sign( ONE, H( I, I+1 ) ) )RESULT( 1+RSUB ) = ULPINV;
            }
         } // 50

         // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)

         // Copy A to VS1, used as workspace

         slacpy(' ', N, N, A, LDA, VS1, LDVS );

         // Compute Q*H and store in HT.

         sgemm('No transpose', 'No transpose', N, N, N, ONE, VS, LDVS, H, LDA, ZERO, HT, LDA );

         // Compute A - Q*H*Q'

         sgemm('No transpose', 'Transpose', N, N, N, -ONE, HT, LDA, VS, LDVS, ONE, VS1, LDVS );

         ANORM = max( SLANGE( '1', N, N, A, LDA, WORK ), SMLNUM );
         WNORM = SLANGE( '1', N, N, VS1, LDVS, WORK );

         if ( ANORM > WNORM ) {
            RESULT[2+RSUB] = ( WNORM / ANORM ) / ( N*ULP );
         } else {
            if ( ANORM < ONE ) {
               RESULT[2+RSUB] = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
            } else {
               RESULT[2+RSUB] = min( WNORM / ANORM, REAL( N ) ) / ( N*ULP );
            }
         }

         // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )

         sort01('Columns', N, N, VS, LDVS, WORK, LWORK, RESULT( 3+RSUB ) );

         // Do Test (4) or Test (10)

         RESULT[4+RSUB] = ZERO;
         for (I = 1; I <= N; I++) { // 60
            if( H( I, I ) != WR( I ) ) RESULT( 4+RSUB ) = ULPINV;
         } // 60
         if ( N > 1 ) {
            if( H( 2, 1 ) == ZERO && WI( 1 ) != ZERO ) RESULT( 4+RSUB ) = ULPINV;
            IF( H( N, N-1 ) == ZERO && WI( N ) != ZERO ) RESULT( 4+RSUB ) = ULPINV;
         }
         for (I = 1; I <= N - 1; I++) { // 70
            if ( H( I+1, I ) != ZERO ) {
               TMP = sqrt( ( H( I+1, I ) ).abs() )* sqrt( ( H( I, I+1 ) ).abs() )                RESULT( 4+RSUB ) = max( RESULT( 4+RSUB ), ABS( WI( I )-TMP ) / max( ULP*TMP, SMLNUM ) )                RESULT( 4+RSUB ) = max( RESULT( 4+RSUB ), ABS( WI( I+1 )+TMP ) / max( ULP*TMP, SMLNUM ) );
            } else if ( I > 1 ) {
               if( H( I+1, I ) == ZERO && H( I, I-1 ) == ZERO && WI( I ) != ZERO )RESULT( 4+RSUB ) = ULPINV;
            }
         } // 70

         // Do Test (5) or Test (11)

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('N', SORT, SSLECT, 'N', N, HT, LDA, SDIM, WRT, WIT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[5+RSUB] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX2', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX2', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         RESULT[5+RSUB] = ZERO;
         for (J = 1; J <= N; J++) { // 90
            for (I = 1; I <= N; I++) { // 80
               if( H( I, J ) != HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV;
            } // 80
         } // 90

         // Do Test (6) or Test (12)

         RESULT[6+RSUB] = ZERO;
         for (I = 1; I <= N; I++) { // 100
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 6+RSUB ) = ULPINV;
         } // 100

         // Do Test (13)

         if ( ISORT == 1 ) {
            RESULT[13] = ZERO;
            KNTEIG = 0;
            for (I = 1; I <= N; I++) { // 110
               if( SSLECT( WR( I ), WI( I ) ) || SSLECT( WR( I ), -WI( I ) ) )KNTEIG = KNTEIG + 1;
               if ( I < N ) {
                  if( ( SSLECT( WR( I+1 ), WI( I+1 ) ) || SSLECT( WR( I+1 ), -WI( I+1 ) ) ) && ( !( SSLECT( WR( I ), WI( I ) ) || SSLECT( WR( I ), -WI( I ) ) ) ) && IINFO != N+2 )RESULT( 13 ) = ULPINV;
               }
            } // 110
            if (SDIM != KNTEIG) RESULT( 13 ) = ULPINV;
         }

      } // 120

      // If there is enough workspace, perform tests (14) and (15)
      // as well as (10) through (13)

      if ( LWORK >= N+( N*N ) / 2 ) {

         // Compute both RCONDE and RCONDV with VS

         SORT = 'S';
         RESULT[14] = ZERO;
         RESULT[15] = ZERO;
         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('V', SORT, SSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX3', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX3', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 140
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 130
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 130
         } // 140
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute both RCONDE and RCONDV without VS, and compare

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('N', SORT, SSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX4', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX4', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Perform tests (14) and (15)

         if (RCNDE1 != RCONDE) RESULT( 14 ) = ULPINV;
         IF( RCNDV1 != RCONDV ) RESULT( 15 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 160
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 150
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 150
         } // 160
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDE with VS, and compare

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('V', SORT, SSLECT, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[14] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX5', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX5', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Perform test (14)

         if (RCNDE1 != RCONDE) RESULT( 14 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 180
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 170
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 170
         } // 180
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDE without VS, and compare

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('N', SORT, SSLECT, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[14] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX6', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX6', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Perform test (14)

         if (RCNDE1 != RCONDE) RESULT( 14 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 200
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 190
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 190
         } // 200
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDV with VS, and compare

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('V', SORT, SSLECT, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX7', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX7', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Perform test (15)

         if (RCNDV1 != RCONDV) RESULT( 15 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 220
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 210
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 210
         } // 220
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDV without VS, and compare

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('N', SORT, SSLECT, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEESX8', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEESX8', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Perform test (15)

         if (RCNDV1 != RCONDV) RESULT( 15 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 240
            if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 230
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 230
         } // 240
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

      }

      } // 250

      // If there are precomputed reciprocal condition numbers, compare
      // computed values with them.

      if ( COMP ) {

         // First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
         // the logical function SSLECT selects the eigenvalues specified
         // by NSLCT and ISLCT.

         SELDIM = N;
         SELOPT = 1;
         EPS = max( ULP, EPSIN );
         for (I = 1; I <= N; I++) { // 260
            IPNT[I] = I;
            SELVAL[I] = false;
            SELWR[I] = WRTMP( I );
            SELWI[I] = WITMP( I );
         } // 260
         for (I = 1; I <= N - 1; I++) { // 280
            KMIN = I;
            VRMIN = WRTMP( I );
            VIMIN = WITMP( I );
            for (J = I + 1; J <= N; J++) { // 270
               if ( WRTMP( J ) < VRMIN ) {
                  KMIN = J;
                  VRMIN = WRTMP( J );
                  VIMIN = WITMP( J );
               }
            } // 270
            WRTMP[KMIN] = WRTMP( I );
            WITMP[KMIN] = WITMP( I );
            WRTMP[I] = VRMIN;
            WITMP[I] = VIMIN;
            ITMP = IPNT( I );
            IPNT[I] = IPNT( KMIN );
            IPNT[KMIN] = ITMP;
         } // 280
         for (I = 1; I <= NSLCT; I++) { // 290
            SELVAL[IPNT( ISLCT( I ) )] = true;
         } // 290

         // Compute condition numbers

         slacpy('F', N, N, A, LDA, HT, LDA );
         sgeesx('N', 'S', SSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO );
         if ( IINFO != 0 && IINFO != N+2 ) {
            RESULT[16] = ULPINV;
            RESULT[17] = ULPINV;
            WRITE( NOUNIT, FMT = 9999 )'SGEESX9', IINFO, N, ISEED( 1 );
            INFO = ( IINFO ).abs();
            GO TO 300;
         }

         // Compare condition number for average of selected eigenvalues
         // taking its condition number into account

         ANORM = SLANGE( '1', N, N, A, LDA, WORK );
         V = max( double( N )*EPS*ANORM, SMLNUM );
         if (ANORM == ZERO) V = ONE;
         if ( V > RCONDV ) {
            TOL = ONE;
         } else {
            TOL = V / RCONDV;
         }
         if ( V > RCDVIN ) {
            TOLIN = ONE;
         } else {
            TOLIN = V / RCDVIN;
         }
         TOL = max( TOL, SMLNUM / EPS );
         TOLIN = max( TOLIN, SMLNUM / EPS );
         if ( EPS*( RCDEIN-TOLIN ) > RCONDE+TOL ) {
            RESULT[16] = ULPINV;
         } else if ( RCDEIN-TOLIN > RCONDE+TOL ) {
            RESULT[16] = ( RCDEIN-TOLIN ) / ( RCONDE+TOL );
         } else if ( RCDEIN+TOLIN < EPS*( RCONDE-TOL ) ) {
            RESULT[16] = ULPINV;
         } else if ( RCDEIN+TOLIN < RCONDE-TOL ) {
            RESULT[16] = ( RCONDE-TOL ) / ( RCDEIN+TOLIN );
         } else {
            RESULT[16] = ONE;
         }

         // Compare condition numbers for right invariant subspace
         // taking its condition number into account

         if ( V > RCONDV*RCONDE ) {
            TOL = RCONDV;
         } else {
            TOL = V / RCONDE;
         }
         if ( V > RCDVIN*RCDEIN ) {
            TOLIN = RCDVIN;
         } else {
            TOLIN = V / RCDEIN;
         }
         TOL = max( TOL, SMLNUM / EPS );
         TOLIN = max( TOLIN, SMLNUM / EPS );
         if ( EPS*( RCDVIN-TOLIN ) > RCONDV+TOL ) {
            RESULT[17] = ULPINV;
         } else if ( RCDVIN-TOLIN > RCONDV+TOL ) {
            RESULT[17] = ( RCDVIN-TOLIN ) / ( RCONDV+TOL );
         } else if ( RCDVIN+TOLIN < EPS*( RCONDV-TOL ) ) {
            RESULT[17] = ULPINV;
         } else if ( RCDVIN+TOLIN < RCONDV-TOL ) {
            RESULT[17] = ( RCONDV-TOL ) / ( RCDVIN+TOLIN );
         } else {
            RESULT[17] = ONE;
         }

         } // 300

      }

 9999 FORMAT( ' SGET24: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, INPUT EXAMPLE NUMBER = ${.i4}');
 9998 FORMAT( ' SGET24: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );

      }
