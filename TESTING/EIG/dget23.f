      SUBROUTINE DGET23( COMP, BALANC, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, LWORK, IWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               COMP;
      String             BALANC;
      int                INFO, JTYPE, LDA, LDLRE, LDVL, LDVR, LWORK, N, NOUNIT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 ), IWORK( * );
      double             A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), SCALE( * ), SCALE1( * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * );
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 ;
      double             EPSIN;
      const              EPSIN = 5.9605D-8 ;
      // ..
      // .. Local Scalars ..
      bool               BALOK, NOBAL;
      String             SENSE;
      int                I, IHI, IHI1, IINFO, ILO, ILO1, ISENS, ISENSM, J, JJ, KMIN;
      double             ABNRM, ABNRM1, EPS, SMLNUM, TNRM, TOL, TOLIN, ULP, ULPINV, V, VIMIN, VMAX, VMX, VRMIN, VRMX, VTST;
      // ..
      // .. Local Arrays ..
      String             SENS( 2 );
      double             DUM( 1 ), RES( 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLAPY2, DNRM2;
      // EXTERNAL LSAME, DLAMCH, DLAPY2, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEEVX, DGET22, DLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Data statements ..
      DATA               SENS / 'N', 'V' /
      // ..
      // .. Executable Statements ..

      // Check for errors

      NOBAL = LSAME( BALANC, 'N' )
      BALOK = NOBAL .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'B' )
      INFO = 0
      if ( .NOT.BALOK ) {
         INFO = -2
      } else if ( THRESH.LT.ZERO ) {
         INFO = -4
      } else if ( NOUNIT.LE.0 ) {
         INFO = -6
      } else if ( N.LT.0 ) {
         INFO = -7
      } else if ( LDA.LT.1 .OR. LDA.LT.N ) {
         INFO = -9
      } else if ( LDVL.LT.1 .OR. LDVL.LT.N ) {
         INFO = -16
      } else if ( LDVR.LT.1 .OR. LDVR.LT.N ) {
         INFO = -18
      } else if ( LDLRE.LT.1 .OR. LDLRE.LT.N ) {
         INFO = -20
      } else if ( LWORK.LT.3*N .OR. ( COMP .AND. LWORK.LT.6*N+N*N ) ) {
         INFO = -31
      }

      if ( INFO.NE.0 ) {
         xerbla('DGET23', -INFO );
         RETURN
      }

      // Quick return if nothing to do

      for (I = 1; I <= 11; I++) { // 10
         RESULT( I ) = -ONE
      } // 10

      IF( N.EQ.0 ) RETURN

      // More Important constants

      ULP = DLAMCH( 'Precision' )
      SMLNUM = DLAMCH( 'S' )
      ULPINV = ONE / ULP

      // Compute eigenvalues and eigenvectors, and test them

      if ( LWORK.GE.6*N+N*N ) {
         SENSE = 'B'
         ISENSM = 2
      } else {
         SENSE = 'E'
         ISENSM = 1
      }
      dlacpy('F', N, N, A, LDA, H, LDA );
      dgeevx(BALANC, 'V', 'V', SENSE, N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, IINFO );
      if ( IINFO.NE.0 ) {
         RESULT( 1 ) = ULPINV
         if ( JTYPE.NE.22 ) {
            WRITE( NOUNIT, FMT = 9998 )'DGEEVX1', IINFO, N, JTYPE, BALANC, ISEED
         } else {
            WRITE( NOUNIT, FMT = 9999 )'DGEEVX1', IINFO, N, ISEED( 1 )
         }
         INFO = ABS( IINFO )
         RETURN
      }

      // Do Test (1)

      dget22('N', 'N', 'N', N, A, LDA, VR, LDVR, WR, WI, WORK, RES );
      RESULT( 1 ) = RES( 1 )

      // Do Test (2)

      dget22('T', 'N', 'T', N, A, LDA, VL, LDVL, WR, WI, WORK, RES );
      RESULT( 2 ) = RES( 1 )

      // Do Test (3)

      for (J = 1; J <= N; J++) { // 30
         TNRM = ONE
         if ( WI( J ).EQ.ZERO ) {
            TNRM = DNRM2( N, VR( 1, J ), 1 )
         } else if ( WI( J ).GT.ZERO ) {
            TNRM = DLAPY2( DNRM2( N, VR( 1, J ), 1 ), DNRM2( N, VR( 1, J+1 ), 1 ) )
         }
         RESULT( 3 ) = MAX( RESULT( 3 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
         if ( WI( J ).GT.ZERO ) {
            VMX = ZERO
            VRMX = ZERO
            for (JJ = 1; JJ <= N; JJ++) { // 20
               VTST = DLAPY2( VR( JJ, J ), VR( JJ, J+1 ) )
               IF( VTST.GT.VMX ) VMX = VTST                IF( VR( JJ, J+1 ).EQ.ZERO .AND. ABS( VR( JJ, J ) ).GT. VRMX )VRMX = ABS( VR( JJ, J ) )
            } // 20
            IF( VRMX / VMX.LT.ONE-TWO*ULP ) RESULT( 3 ) = ULPINV
         }
      } // 30

      // Do Test (4)

      for (J = 1; J <= N; J++) { // 50
         TNRM = ONE
         if ( WI( J ).EQ.ZERO ) {
            TNRM = DNRM2( N, VL( 1, J ), 1 )
         } else if ( WI( J ).GT.ZERO ) {
            TNRM = DLAPY2( DNRM2( N, VL( 1, J ), 1 ), DNRM2( N, VL( 1, J+1 ), 1 ) )
         }
         RESULT( 4 ) = MAX( RESULT( 4 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
         if ( WI( J ).GT.ZERO ) {
            VMX = ZERO
            VRMX = ZERO
            for (JJ = 1; JJ <= N; JJ++) { // 40
               VTST = DLAPY2( VL( JJ, J ), VL( JJ, J+1 ) )
               IF( VTST.GT.VMX ) VMX = VTST                IF( VL( JJ, J+1 ).EQ.ZERO .AND. ABS( VL( JJ, J ) ).GT. VRMX )VRMX = ABS( VL( JJ, J ) )
            } // 40
            IF( VRMX / VMX.LT.ONE-TWO*ULP ) RESULT( 4 ) = ULPINV
         }
      } // 50

      // Test for all options of computing condition numbers

      for (ISENS = 1; ISENS <= ISENSM; ISENS++) { // 200

         SENSE = SENS( ISENS )

         // Compute eigenvalues only, and test them

         dlacpy('F', N, N, A, LDA, H, LDA );
         dgeevx(BALANC, 'N', 'N', SENSE, N, H, LDA, WR1, WI1, DUM, 1, DUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, IWORK, IINFO );
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEEVX2', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEEVX2', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5)

         for (J = 1; J <= N; J++) { // 60
            IF( WR( J ).NE.WR1( J ) .OR. WI( J ).NE.WI1( J ) ) RESULT( 5 ) = ULPINV
         } // 60

         // Do Test (8)

         if ( .NOT.NOBAL ) {
            for (J = 1; J <= N; J++) { // 70
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
            } // 70
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         }

         // Do Test (9)

         if ( ISENS.EQ.2 .AND. N.GT.1 ) {
            for (J = 1; J <= N; J++) { // 80
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
            } // 80
         }

         // Compute eigenvalues and right eigenvectors, and test them

         dlacpy('F', N, N, A, LDA, H, LDA );
         dgeevx(BALANC, 'N', 'V', SENSE, N, H, LDA, WR1, WI1, DUM, 1, LRE, LDLRE, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, IWORK, IINFO );
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEEVX3', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEEVX3', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 90
            IF( WR( J ).NE.WR1( J ) .OR. WI( J ).NE.WI1( J ) ) RESULT( 5 ) = ULPINV
         } // 90

         // Do Test (6)

         for (J = 1; J <= N; J++) { // 110
            for (JJ = 1; JJ <= N; JJ++) { // 100
               IF( VR( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 6 ) = ULPINV
            } // 100
         } // 110

         // Do Test (8) again

         if ( .NOT.NOBAL ) {
            for (J = 1; J <= N; J++) { // 120
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
            } // 120
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         }

         // Do Test (9) again

         if ( ISENS.EQ.2 .AND. N.GT.1 ) {
            for (J = 1; J <= N; J++) { // 130
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
            } // 130
         }

         // Compute eigenvalues and left eigenvectors, and test them

         dlacpy('F', N, N, A, LDA, H, LDA );
         dgeevx(BALANC, 'V', 'N', SENSE, N, H, LDA, WR1, WI1, LRE, LDLRE, DUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, IWORK, IINFO );
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEEVX4', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEEVX4', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 140
            IF( WR( J ).NE.WR1( J ) .OR. WI( J ).NE.WI1( J ) ) RESULT( 5 ) = ULPINV
         } // 140

         // Do Test (7)

         for (J = 1; J <= N; J++) { // 160
            for (JJ = 1; JJ <= N; JJ++) { // 150
               IF( VL( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 7 ) = ULPINV
            } // 150
         } // 160

         // Do Test (8) again

         if ( .NOT.NOBAL ) {
            for (J = 1; J <= N; J++) { // 170
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
            } // 170
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         }

         // Do Test (9) again

         if ( ISENS.EQ.2 .AND. N.GT.1 ) {
            for (J = 1; J <= N; J++) { // 180
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
            } // 180
         }

         } // 190

      } // 200

      // If COMP, compare condition numbers to precomputed ones

      if ( COMP ) {
         dlacpy('F', N, N, A, LDA, H, LDA );
         dgeevx('N', 'V', 'V', 'B', N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, IINFO );
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'DGEEVX5', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         DO 220 I = 1, N - 1
            KMIN = I
            VRMIN = WR( I )
            VIMIN = WI( I )
            DO 210 J = I + 1, N
               if ( WR( J ).LT.VRMIN ) {
                  KMIN = J
                  VRMIN = WR( J )
                  VIMIN = WI( J )
               }
            } // 210
            WR( KMIN ) = WR( I )
            WI( KMIN ) = WI( I )
            WR( I ) = VRMIN
            WI( I ) = VIMIN
            VRMIN = RCONDE( KMIN )
            RCONDE( KMIN ) = RCONDE( I )
            RCONDE( I ) = VRMIN
            VRMIN = RCONDV( KMIN )
            RCONDV( KMIN ) = RCONDV( I )
            RCONDV( I ) = VRMIN
         } // 220

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         RESULT( 10 ) = ZERO
         EPS = MAX( EPSIN, ULP )
         V = MAX( DBLE( N )*EPS*ABNRM, SMLNUM )
         IF( ABNRM.EQ.ZERO ) V = ONE
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
            } else if ( RCDVIN( I )+TOLIN.LT.EPS*( RCONDV( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( RCDVIN( I )+TOLIN.LT.RCONDV( I )-TOL ) {
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
            } else if ( RCDEIN( I )+TOLIN.LT.EPS*( RCONDE( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( RCDEIN( I )+TOLIN.LT.RCONDE( I )-TOL ) {
               VMAX = ( RCONDE( I )-TOL ) / ( RCDEIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            RESULT( 11 ) = MAX( RESULT( 11 ), VMAX )
         } // 240
         } // 250

      }

 9999 FORMAT( ' DGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' DGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', BALANC = ', A, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of DGET23

      }
