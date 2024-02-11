import 'package:lapack/src/complex.dart';

import 'common.dart';

      void zget24(final int COMP, final int JTYPE, final int THRESH, final Array<int> ISEED, final int NOUNIT, final int N, final Matrix<double> A, final int LDA, final int H, final int HT, final int W, final int WT, final int WTMP, final Matrix<double> VS, final int LDVS, final int VS1, final int RCDEIN, final int RCDVIN, final int NSLCT, final int ISLCT, final int ISRT, final int RESULT, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Array<bool> BWORK, final Box<int> INFO,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               COMP;
      int                INFO, ISRT, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT;
      double             RCDEIN, RCDVIN, THRESH;
      bool               BWORK( * );
      int                ISEED( 4 ), ISLCT( * );
      double             RESULT( 17 ), RWORK( * );
      Complex         A( LDA, * ), H( LDA, * ), HT( LDA, * ), VS( LDVS, * ), VS1( LDVS, * ), W( * ), WORK( * ), WT( * ), WTMP( * );
      // ..

      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             EPSIN;
      const              EPSIN = 5.9605e-8 ;
      String             SORT;
      int                I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, RSUB, SDIM, SDIM1;
      double             ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, SMLNUM, TOL, TOLIN, ULP, ULPINV, V, VRICMP, VRIMIN, WNORM;
      Complex         CTMP;
      int                IPNT( 20 );
      // ..
      // .. External Functions ..
      //- bool               ZSLECT;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL ZSLECT, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGEESX, ZGEMM, ZLACPY, ZUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Arrays in Common ..
      // bool               sslct.SELVAL( 20 );
      // double             sslct.SELWI( 20 ), sslct.SELWR( 20 );
      // // ..
      // // .. Scalars in Common ..
      // int                sslct.SELDIM, sslct.SELOPT;
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / sslct.SELOPT, sslct.SELDIM, sslct.SELVAL, sslct.SELWR, sslct.SELWI

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
         xerbla('ZGET24', -INFO );
         return;
      }

      // Quick return if nothing to do

      for (I = 1; I <= 17; I++) { // 10
         RESULT[I] = -ONE;
      } // 10

      if (N == 0) return;

      // Important constants

      SMLNUM = dlamch( 'Safe minimum' );
      ULP = dlamch( 'Precision' );
      ULPINV = ONE / ULP;

      // Perform tests (1)-(13)

      sslct.SELOPT = 0;
      for (ISORT = 0; ISORT <= 1; ISORT++) { // 90
         if ( ISORT == 0 ) {
            SORT = 'N';
            RSUB = 0;
         } else {
            SORT = 'S';
            RSUB = 6;
         }

         // Compute Schur form and Schur vectors, and test them

         zlacpy('F', N, N, A, LDA, H, LDA );
         zgeesx('V', SORT, ZSLECT, 'N', N, H, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1+RSUB] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX1', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX1', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            return;
         }
         if ( ISORT == 0 ) {
            zcopy(N, W, 1, WTMP, 1 );
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

         zlacpy(' ', N, N, A, LDA, VS1, LDVS );

         // Compute Q*H and store in HT.

         zgemm('No transpose', 'No transpose', N, N, N, CONE, VS, LDVS, H, LDA, CZERO, HT, LDA );

         // Compute A - Q*H*Q'

         zgemm('No transpose', 'Conjugate transpose', N, N, N, -CONE, HT, LDA, VS, LDVS, CONE, VS1, LDVS );

         ANORM = max( ZLANGE( '1', N, N, A, LDA, RWORK ), SMLNUM );
         WNORM = ZLANGE( '1', N, N, VS1, LDVS, RWORK );

         if ( ANORM > WNORM ) {
            RESULT[2+RSUB] = ( WNORM / ANORM ) / ( N*ULP );
         } else {
            if ( ANORM < ONE ) {
               RESULT[2+RSUB] = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
            } else {
               RESULT[2+RSUB] = min( WNORM / ANORM, N.toDouble() ) / ( N*ULP );
            }
         }

         // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )

         zunt01('Columns', N, N, VS, LDVS, WORK, LWORK, RWORK, RESULT( 3+RSUB ) );

         // Do Test (4) or Test (10)

         RESULT[4+RSUB] = ZERO;
         for (I = 1; I <= N; I++) { // 40
            if( H( I, I ) != W( I ) ) RESULT( 4+RSUB ) = ULPINV;
         } // 40

         // Do Test (5) or Test (11)

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('N', SORT, ZSLECT, 'N', N, HT, LDA, SDIM, WT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[5+RSUB] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX2', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX2', IINFO, N, ISEED( 1 );
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
               if( ZSLECT( W( I ) ) ) KNTEIG = KNTEIG + 1;
               if ( I < N ) {
                  if[ZSLECT( W( I+1 ) ) && ( !ZSLECT( W( I ) ) ) )RESULT( 13] = ULPINV;
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
         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('V', SORT, ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX3', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX3', IINFO, N, ISEED( 1 );
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

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('N', SORT, ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX4', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX4', IINFO, N, ISEED( 1 );
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

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('V', SORT, ZSLECT, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX5', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX5', IINFO, N, ISEED( 1 );
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

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('N', SORT, ZSLECT, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[14] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX6', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX6', IINFO, N, ISEED( 1 );
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

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('V', SORT, ZSLECT, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX7', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX7', IINFO, N, ISEED( 1 );
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

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('N', SORT, ZSLECT, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[15] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX8', IINFO, N, JTYPE, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX8', IINFO, N, ISEED( 1 );
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

         // First set up sslct.SELOPT, sslct.SELDIM, sslct.SELVAL, sslct.SELWR and sslct.SELWI so that
         // the logical function ZSLECT selects the eigenvalues specified
         // by NSLCT, ISLCT and ISRT.

         sslct.SELDIM = N;
         sslct.SELOPT = 1;
         EPS = max( ULP, EPSIN );
         for (I = 1; I <= N; I++) { // 230
            IPNT[I] = I;
            sslct.SELVAL[I] = false;
            sslct.SELWR[I] = (WTMP( I )).toDouble();
            sslct.SELWI[I] = DIMAG( WTMP( I ) );
         } // 230
         for (I = 1; I <= N - 1; I++) { // 250
            KMIN = I;
            if ( ISRT == 0 ) {
               VRIMIN = (WTMP( I )).toDouble();
            } else {
               VRIMIN = DIMAG( WTMP( I ) );
            }
            for (J = I + 1; J <= N; J++) { // 240
               if ( ISRT == 0 ) {
                  VRICMP = (WTMP( J )).toDouble();
               } else {
                  VRICMP = DIMAG( WTMP( J ) );
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
            sslct.SELVAL[IPNT( ISLCT( I ) )] = true;
         } // 260

         // Compute condition numbers

         zlacpy('F', N, N, A, LDA, HT, LDA );
         zgeesx('N', 'S', ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[16] = ULPINV;
            RESULT[17] = ULPINV;
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX9', IINFO, N, ISEED( 1 );
            INFO = ( IINFO ).abs();
            GO TO 270;
         }

         // Compare condition number for average of selected eigenvalues
         // taking its condition number into account

         ANORM = ZLANGE( '1', N, N, A, LDA, RWORK );
         V = max( N.toDouble()*EPS*ANORM, SMLNUM );
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

 9999 FORMAT( ' ZGET24: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, INPUT EXAMPLE NUMBER = ${.i4}');
 9998 FORMAT( ' ZGET24: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );

      }
