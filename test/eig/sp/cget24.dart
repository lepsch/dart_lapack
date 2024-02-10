      void cget24(COMP, JTYPE, THRESH, final Array<int> ISEED, NOUNIT, N, final Matrix<double> A, final int LDA, H, HT, W, WT, WTMP, final Matrix<double> VS, final int LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, ISRT, RESULT, final Array<double> WORK, final int LWORK, RWORK, BWORK, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               COMP;
      int                INFO, ISRT, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT;
      double               RCDEIN, RCDVIN, THRESH;
      bool               BWORK( * );
      int                ISEED( 4 ), ISLCT( * );
      double               RESULT( 17 ), RWORK( * );
      Complex            A( LDA, * ), H( LDA, * ), HT( LDA, * ), VS( LDVS, * ), VS1( LDVS, * ), W( * ), WORK( * ), WT( * ), WTMP( * );
      // ..

      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               EPSIN;
      const              EPSIN = 5.9605e-8 ;
      String             SORT;
      int                I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, RSUB, SDIM, SDIM1;
      double               ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, SMLNUM, TOL, TOLIN, ULP, ULPINV, V, VRICMP, VRIMIN, WNORM;
      Complex            CTMP;
      int                IPNT( 20 );
      // ..
      // .. External Functions ..
      //- bool               CSLECT;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CSLECT, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEESX, CGEMM, CLACPY, CUNT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
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
         INFO = -15;
      } else if ( LWORK < 2*N ) {
         INFO = -24;
      }

      if ( INFO != 0 ) {
         xerbla('CGET24', -INFO );
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
      for (ISORT = 0; ISORT <= 1; ISORT++) { // 90
         if ( ISORT == 0 ) {
            SORT = 'N';
            RSUB = 0;
         } else {
            SORT = 'S';
            RSUB = 6;
         }

         // Compute Schur form and Schur vectors, and test them

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeesx('V', SORT, CSLECT, 'N', N, H, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1+RSUB] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX1', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX1', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            return;
         }
         if ( ISORT == 0 ) {
            ccopy(N, W, 1, WTMP, 1 );
         }

         // Do Test (1) or Test (7)

         RESULT[1+RSUB] = ZERO;
         for (J = 1; J <= N - 1; J++) { // 30
            for (I = J + 1; I <= N; I++) { // 20
               if( H( I, J ) != CZERO ) RESULT( 1+RSUB ) = ULPINV;
            } // 20
         } // 30

         // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)

         // Copy A to VS1, used as workspace

         clacpy(' ', N, N, A, LDA, VS1, LDVS );

         // Compute Q*H and store in HT.

         cgemm('No transpose', 'No transpose', N, N, N, CONE, VS, LDVS, H, LDA, CZERO, HT, LDA );

         // Compute A - Q*H*Q'

         cgemm('No transpose', 'Conjugate transpose', N, N, N, -CONE, HT, LDA, VS, LDVS, CONE, VS1, LDVS );

         ANORM = max( CLANGE( '1', N, N, A, LDA, RWORK ), SMLNUM );
         WNORM = CLANGE( '1', N, N, VS1, LDVS, RWORK );

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

         cunt01('Columns', N, N, VS, LDVS, WORK, LWORK, RWORK, RESULT( 3+RSUB ) );

         // Do Test (4) or Test (10)

         RESULT[4+RSUB] = ZERO;
         for (I = 1; I <= N; I++) { // 40
            if( H( I, I ) != W( I ) ) RESULT( 4+RSUB ) = ULPINV;
         } // 40

         // Do Test (5) or Test (11)

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('N', SORT, CSLECT, 'N', N, HT, LDA, SDIM, WT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[5+RSUB] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX2', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX2', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         RESULT[5+RSUB] = ZERO;
         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               if( H( I, J ) != HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV;
            } // 50
         } // 60

         // Do Test (6) or Test (12)

         RESULT[6+RSUB] = ZERO;
         for (I = 1; I <= N; I++) { // 70
            if( W( I ) != WT( I ) ) RESULT( 6+RSUB ) = ULPINV;
         } // 70

         // Do Test (13)

         if ( ISORT == 1 ) {
            RESULT[13] = ZERO;
            KNTEIG = 0;
            for (I = 1; I <= N; I++) { // 80
               if( CSLECT( W( I ) ) ) KNTEIG = KNTEIG + 1;
               if ( I < N ) {
                  if[CSLECT( W( I+1 ) ) && ( !CSLECT( W( I ) ) ) )RESULT( 13] = ULPINV;
               }
            } // 80
            if (SDIM != KNTEIG) RESULT( 13 ) = ULPINV;
         }

      } // 90

      // If there is enough workspace, perform tests (14) and (15)
      // as well as (10) through (13)

      if ( LWORK >= ( N*( N+1 ) ) / 2 ) {

         // Compute both RCONDE and RCONDV with VS

         SORT = 'S';
         RESULT[14] = ZERO;
         RESULT[15] = ZERO;
         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('V', SORT, CSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX3', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX3', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 110
            if( W( I ) != WT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 100
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 100
         } // 110
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute both RCONDE and RCONDV without VS, and compare

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('N', SORT, CSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX4', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX4', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         // Perform tests (14) and (15)

         if (RCNDE1 != RCONDE) RESULT( 14 ) = ULPINV;
         IF( RCNDV1 != RCONDV ) RESULT( 15 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 130
            if( W( I ) != WT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 120
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 120
         } // 130
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDE with VS, and compare

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('V', SORT, CSLECT, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX5', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX5', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         // Perform test (14)

         if (RCNDE1 != RCONDE) RESULT( 14 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 150
            if( W( I ) != WT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 140
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 140
         } // 150
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDE without VS, and compare

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('N', SORT, CSLECT, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX6', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX6', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         // Perform test (14)

         if (RCNDE1 != RCONDE) RESULT( 14 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 170
            if( W( I ) != WT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 160
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 160
         } // 170
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDV with VS, and compare

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('V', SORT, CSLECT, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX7', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX7', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         // Perform test (15)

         if (RCNDV1 != RCONDV) RESULT( 15 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 190
            if( W( I ) != WT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 180
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 180
         } // 190
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

         // Compute RCONDV without VS, and compare

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('N', SORT, CSLECT, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEESX8', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEESX8', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 220;
         }

         // Perform test (15)

         if (RCNDV1 != RCONDV) RESULT( 15 ) = ULPINV;

         // Perform tests (10), (11), (12), and (13)

         for (I = 1; I <= N; I++) { // 210
            if( W( I ) != WT( I ) ) RESULT( 10 ) = ULPINV;
            for (J = 1; J <= N; J++) { // 200
               if( H( I, J ) != HT( I, J ) ) RESULT( 11 ) = ULPINV;
               IF( VS( I, J ) != VS1( I, J ) ) RESULT( 12 ) = ULPINV;
            } // 200
         } // 210
         if (SDIM != SDIM1) RESULT( 13 ) = ULPINV;

      }

      } // 220

      // If there are precomputed reciprocal condition numbers, compare
      // computed values with them.

      if ( COMP ) {

         // First set up SELOPT, SELDIM, SELVAL, SELWR and SELWI so that
         // the logical function CSLECT selects the eigenvalues specified
         // by NSLCT, ISLCT and ISRT.

         SELDIM = N;
         SELOPT = 1;
         EPS = max( ULP, EPSIN );
         for (I = 1; I <= N; I++) { // 230
            IPNT[I] = I;
            SELVAL[I] = false;
            SELWR[I] = double( WTMP( I ) );
            SELWI[I] = AIMAG( WTMP( I ) );
         } // 230
         for (I = 1; I <= N - 1; I++) { // 250
            KMIN = I;
            if ( ISRT == 0 ) {
               VRIMIN = double( WTMP( I ) );
            } else {
               VRIMIN = AIMAG( WTMP( I ) );
            }
            for (J = I + 1; J <= N; J++) { // 240
               if ( ISRT == 0 ) {
                  VRICMP = double( WTMP( J ) );
               } else {
                  VRICMP = AIMAG( WTMP( J ) );
               }
               if ( VRICMP < VRIMIN ) {
                  KMIN = J;
                  VRIMIN = VRICMP;
               }
            } // 240
            CTMP = WTMP( KMIN );
            WTMP[KMIN] = WTMP( I );
            WTMP[I] = CTMP;
            ITMP = IPNT( I );
            IPNT[I] = IPNT( KMIN );
            IPNT[KMIN] = ITMP;
         } // 250
         for (I = 1; I <= NSLCT; I++) { // 260
            SELVAL[IPNT( ISLCT( I ) )] = true;
         } // 260

         // Compute condition numbers

         clacpy('F', N, N, A, LDA, HT, LDA );
         cgeesx('N', 'S', CSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[16] = ULPINV;
            RESULT[17] = ULPINV;
            WRITE( NOUNIT, FMT = 9999 )'CGEESX9', IINFO, N, ISEED( 1 );
            INFO = ( IINFO ).abs();
            GO TO 270;
         }

         // Compare condition number for average of selected eigenvalues
         // taking its condition number into account

         ANORM = CLANGE( '1', N, N, A, LDA, RWORK );
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

         } // 270

      }

 9999 FORMAT( ' CGET24: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, INPUT EXAMPLE NUMBER = ${.i4}');
 9998 FORMAT( ' CGET24: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );

      }
