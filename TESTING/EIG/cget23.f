      SUBROUTINE CGET23( COMP, ISRT, BALANC, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, H, W, W1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, LWORK, RWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               COMP;
      String             BALANC;
      int                INFO, ISRT, JTYPE, LDA, LDLRE, LDVL, LDVR, LWORK, N, NOUNIT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), RWORK( * ), SCALE( * ), SCALE1( * )
      COMPLEX            A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), W1( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 ;
      REAL               EPSIN
      const              EPSIN = 5.9605E-8 ;
      // ..
      // .. Local Scalars ..
      bool               BALOK, NOBAL;
      String             SENSE;
      int                I, IHI, IHI1, IINFO, ILO, ILO1, ISENS, ISENSM, J, JJ, KMIN;
      REAL               ABNRM, ABNRM1, EPS, SMLNUM, TNRM, TOL, TOLIN, ULP, ULPINV, V, VMAX, VMX, VRICMP, VRIMIN, VRMX, VTST;
      COMPLEX            CTMP
      // ..
      // .. Local Arrays ..
      String             SENS( 2 );
      REAL               RES( 2 )
      COMPLEX            CDUM( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SCNRM2, SLAMCH
      // EXTERNAL LSAME, SCNRM2, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEEVX, CGET22, CLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Data statements ..
      DATA               SENS / 'N', 'V' /
      // ..
      // .. Executable Statements ..

      // Check for errors

      NOBAL = LSAME( BALANC, 'N' )
      BALOK = NOBAL || LSAME( BALANC, 'P' ) || LSAME( BALANC, 'S' ) || LSAME( BALANC, 'B' )
      INFO = 0
      if ( ISRT != 0 && ISRT != 1 ) {
         INFO = -2
      } else if ( .NOT.BALOK ) {
         INFO = -3
      } else if ( THRESH < ZERO ) {
         INFO = -5
      } else if ( NOUNIT.LE.0 ) {
         INFO = -7
      } else if ( N < 0 ) {
         INFO = -8
      } else if ( LDA < 1 || LDA < N ) {
         INFO = -10
      } else if ( LDVL < 1 || LDVL < N ) {
         INFO = -15
      } else if ( LDVR < 1 || LDVR < N ) {
         INFO = -17
      } else if ( LDLRE < 1 || LDLRE < N ) {
         INFO = -19
      } else if ( LWORK < 2*N || ( COMP && LWORK < 2*N+N*N ) ) {
         INFO = -30
      }

      if ( INFO != 0 ) {
         xerbla('CGET23', -INFO );
         RETURN
      }

      // Quick return if nothing to do

      for (I = 1; I <= 11; I++) { // 10
         RESULT( I ) = -ONE
      } // 10

      if (N == 0) RETURN;

      // More Important constants

      ULP = SLAMCH( 'Precision' )
      SMLNUM = SLAMCH( 'S' )
      ULPINV = ONE / ULP

      // Compute eigenvalues and eigenvectors, and test them

      if ( LWORK.GE.2*N+N*N ) {
         SENSE = 'B'
         ISENSM = 2
      } else {
         SENSE = 'E'
         ISENSM = 1
      }
      clacpy('F', N, N, A, LDA, H, LDA );
      cgeevx(BALANC, 'V', 'V', SENSE, N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO );
      if ( IINFO != 0 ) {
         RESULT( 1 ) = ULPINV
         if ( JTYPE != 22 ) {
            WRITE( NOUNIT, FMT = 9998 )'CGEEVX1', IINFO, N, JTYPE, BALANC, ISEED
         } else {
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX1', IINFO, N, ISEED( 1 )
         }
         INFO = ABS( IINFO )
         RETURN
      }

      // Do Test (1)

      cget22('N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES );
      RESULT( 1 ) = RES( 1 )

      // Do Test (2)

      cget22('C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES );
      RESULT( 2 ) = RES( 1 )

      // Do Test (3)

      for (J = 1; J <= N; J++) { // 30
         TNRM = SCNRM2( N, VR( 1, J ), 1 )
         RESULT( 3 ) = MAX( RESULT( 3 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
         VMX = ZERO
         VRMX = ZERO
         for (JJ = 1; JJ <= N; JJ++) { // 20
            VTST = ABS( VR( JJ, J ) )
            if (VTST.GT.VMX) VMX = VTST             IF( AIMAG( VR( JJ, J ) ) == ZERO && ABS( REAL( VR( JJ, J ) ) ).GT.VRMX ) VRMX = ABS( REAL( VR( JJ, J ) ) );
         } // 20
         if (VRMX / VMX < ONE-TWO*ULP) RESULT( 3 ) = ULPINV;
      } // 30

      // Do Test (4)

      for (J = 1; J <= N; J++) { // 50
         TNRM = SCNRM2( N, VL( 1, J ), 1 )
         RESULT( 4 ) = MAX( RESULT( 4 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
         VMX = ZERO
         VRMX = ZERO
         for (JJ = 1; JJ <= N; JJ++) { // 40
            VTST = ABS( VL( JJ, J ) )
            if (VTST.GT.VMX) VMX = VTST             IF( AIMAG( VL( JJ, J ) ) == ZERO && ABS( REAL( VL( JJ, J ) ) ).GT.VRMX ) VRMX = ABS( REAL( VL( JJ, J ) ) );
         } // 40
         if (VRMX / VMX < ONE-TWO*ULP) RESULT( 4 ) = ULPINV;
      } // 50

      // Test for all options of computing condition numbers

      for (ISENS = 1; ISENS <= ISENSM; ISENS++) { // 200

         SENSE = SENS( ISENS )

         // Compute eigenvalues only, and test them

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx(BALANC, 'N', 'N', SENSE, N, H, LDA, W1, CDUM, 1, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX2', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX2', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5)

         for (J = 1; J <= N; J++) { // 60
            IF( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV
         } // 60

         // Do Test (8)

         if ( .NOT.NOBAL ) {
            for (J = 1; J <= N; J++) { // 70
               IF( SCALE( J ) != SCALE1( J ) ) RESULT( 8 ) = ULPINV
            } // 70
            if (ILO != ILO1) RESULT( 8 ) = ULPINV             IF( IHI != IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM != ABNRM1 ) RESULT( 8 ) = ULPINV;
         }

         // Do Test (9)

         if ( ISENS == 2 && N.GT.1 ) {
            for (J = 1; J <= N; J++) { // 80
               IF( RCONDV( J ) != RCNDV1( J ) ) RESULT( 9 ) = ULPINV
            } // 80
         }

         // Compute eigenvalues and right eigenvectors, and test them

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx(BALANC, 'N', 'V', SENSE, N, H, LDA, W1, CDUM, 1, LRE, LDLRE, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX3', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX3', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 90
            IF( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV
         } // 90

         // Do Test (6)

         for (J = 1; J <= N; J++) { // 110
            for (JJ = 1; JJ <= N; JJ++) { // 100
               IF( VR( J, JJ ) != LRE( J, JJ ) ) RESULT( 6 ) = ULPINV
            } // 100
         } // 110

         // Do Test (8) again

         if ( .NOT.NOBAL ) {
            for (J = 1; J <= N; J++) { // 120
               IF( SCALE( J ) != SCALE1( J ) ) RESULT( 8 ) = ULPINV
            } // 120
            if (ILO != ILO1) RESULT( 8 ) = ULPINV             IF( IHI != IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM != ABNRM1 ) RESULT( 8 ) = ULPINV;
         }

         // Do Test (9) again

         if ( ISENS == 2 && N.GT.1 ) {
            for (J = 1; J <= N; J++) { // 130
               IF( RCONDV( J ) != RCNDV1( J ) ) RESULT( 9 ) = ULPINV
            } // 130
         }

         // Compute eigenvalues and left eigenvectors, and test them

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx(BALANC, 'V', 'N', SENSE, N, H, LDA, W1, LRE, LDLRE, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX4', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX4', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 140
            IF( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV
         } // 140

         // Do Test (7)

         for (J = 1; J <= N; J++) { // 160
            for (JJ = 1; JJ <= N; JJ++) { // 150
               IF( VL( J, JJ ) != LRE( J, JJ ) ) RESULT( 7 ) = ULPINV
            } // 150
         } // 160

         // Do Test (8) again

         if ( .NOT.NOBAL ) {
            for (J = 1; J <= N; J++) { // 170
               IF( SCALE( J ) != SCALE1( J ) ) RESULT( 8 ) = ULPINV
            } // 170
            if (ILO != ILO1) RESULT( 8 ) = ULPINV             IF( IHI != IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM != ABNRM1 ) RESULT( 8 ) = ULPINV;
         }

         // Do Test (9) again

         if ( ISENS == 2 && N.GT.1 ) {
            for (J = 1; J <= N; J++) { // 180
               IF( RCONDV( J ) != RCNDV1( J ) ) RESULT( 9 ) = ULPINV
            } // 180
         }

         } // 190

      } // 200

      // If COMP, compare condition numbers to precomputed ones

      if ( COMP ) {
         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx('N', 'V', 'V', 'B', N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT( 1 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX5', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         for (I = 1; I <= N - 1; I++) { // 220
            KMIN = I
            if ( ISRT == 0 ) {
               VRIMIN = REAL( W( I ) )
            } else {
               VRIMIN = AIMAG( W( I ) )
            }
            for (J = I + 1; J <= N; J++) { // 210
               if ( ISRT == 0 ) {
                  VRICMP = REAL( W( J ) )
               } else {
                  VRICMP = AIMAG( W( J ) )
               }
               if ( VRICMP < VRIMIN ) {
                  KMIN = J
                  VRIMIN = VRICMP
               }
            } // 210
            CTMP = W( KMIN )
            W( KMIN ) = W( I )
            W( I ) = CTMP
            VRIMIN = RCONDE( KMIN )
            RCONDE( KMIN ) = RCONDE( I )
            RCONDE( I ) = VRIMIN
            VRIMIN = RCONDV( KMIN )
            RCONDV( KMIN ) = RCONDV( I )
            RCONDV( I ) = VRIMIN
         } // 220

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         RESULT( 10 ) = ZERO
         EPS = MAX( EPSIN, ULP )
         V = MAX( REAL( N )*EPS*ABNRM, SMLNUM )
         if (ABNRM == ZERO) V = ONE;
         for (I = 1; I <= N; I++) { // 230
            if ( V.GT.RCONDV( I )*RCONDE( I ) ) {
               TOL = RCONDV( I )
            } else {
               TOL = V / RCONDE( I )
            }
            if ( V.GT.RCDVIN( I )*RCDEIN( I ) ) {
               TOLIN = RCDVIN( I )
            } else {
               TOLIN = V / RCDEIN( I )
            }
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            if ( EPS*( RCDVIN( I )-TOLIN ).GT.RCONDV( I )+TOL ) {
               VMAX = ONE / EPS
            } else if ( RCDVIN( I )-TOLIN.GT.RCONDV( I )+TOL ) {
               VMAX = ( RCDVIN( I )-TOLIN ) / ( RCONDV( I )+TOL )
            } else if ( RCDVIN( I )+TOLIN < EPS*( RCONDV( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( RCDVIN( I )+TOLIN < RCONDV( I )-TOL ) {
               VMAX = ( RCONDV( I )-TOL ) / ( RCDVIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            RESULT( 10 ) = MAX( RESULT( 10 ), VMAX )
         } // 230

         // Compare condition numbers for eigenvalues
         // taking their condition numbers into account

         RESULT( 11 ) = ZERO
         for (I = 1; I <= N; I++) { // 240
            if ( V.GT.RCONDV( I ) ) {
               TOL = ONE
            } else {
               TOL = V / RCONDV( I )
            }
            if ( V.GT.RCDVIN( I ) ) {
               TOLIN = ONE
            } else {
               TOLIN = V / RCDVIN( I )
            }
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            if ( EPS*( RCDEIN( I )-TOLIN ).GT.RCONDE( I )+TOL ) {
               VMAX = ONE / EPS
            } else if ( RCDEIN( I )-TOLIN.GT.RCONDE( I )+TOL ) {
               VMAX = ( RCDEIN( I )-TOLIN ) / ( RCONDE( I )+TOL )
            } else if ( RCDEIN( I )+TOLIN < EPS*( RCONDE( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( RCDEIN( I )+TOLIN < RCONDE( I )-TOL ) {
               VMAX = ( RCONDE( I )-TOL ) / ( RCDEIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            RESULT( 11 ) = MAX( RESULT( 11 ), VMAX )
         } // 240
         } // 250

      }

 9999 FORMAT( ' CGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' CGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', BALANC = ', A, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of CGET23

      }
