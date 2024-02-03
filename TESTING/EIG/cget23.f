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
      BALOK = NOBAL .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'B' )
      INFO = 0
      if ( ISRT.NE.0 .AND. ISRT.NE.1 ) {
         INFO = -2
      } else if ( .NOT.BALOK ) {
         INFO = -3
      } else if ( THRESH.LT.ZERO ) {
         INFO = -5
      } else if ( NOUNIT.LE.0 ) {
         INFO = -7
      } else if ( N.LT.0 ) {
         INFO = -8
      } else if ( LDA.LT.1 .OR. LDA.LT.N ) {
         INFO = -10
      } else if ( LDVL.LT.1 .OR. LDVL.LT.N ) {
         INFO = -15
      } else if ( LDVR.LT.1 .OR. LDVR.LT.N ) {
         INFO = -17
      } else if ( LDLRE.LT.1 .OR. LDLRE.LT.N ) {
         INFO = -19
      } else if ( LWORK.LT.2*N .OR. ( COMP .AND. LWORK.LT.2*N+N*N ) ) {
         INFO = -30
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGET23', -INFO )
         RETURN
      }

      // Quick return if nothing to do

      DO 10 I = 1, 11
         RESULT( I ) = -ONE
   10 CONTINUE

      IF( N.EQ.0 ) RETURN

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
      CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
      CALL CGEEVX( BALANC, 'V', 'V', SENSE, N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO )
      if ( IINFO.NE.0 ) {
         RESULT( 1 ) = ULPINV
         if ( JTYPE.NE.22 ) {
            WRITE( NOUNIT, FMT = 9998 )'CGEEVX1', IINFO, N, JTYPE, BALANC, ISEED
         } else {
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX1', IINFO, N, ISEED( 1 )
         }
         INFO = ABS( IINFO )
         RETURN
      }

      // Do Test (1)

      CALL CGET22( 'N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES )
      RESULT( 1 ) = RES( 1 )

      // Do Test (2)

      CALL CGET22( 'C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES )
      RESULT( 2 ) = RES( 1 )

      // Do Test (3)

      DO 30 J = 1, N
         TNRM = SCNRM2( N, VR( 1, J ), 1 )
         RESULT( 3 ) = MAX( RESULT( 3 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
         VMX = ZERO
         VRMX = ZERO
         DO 20 JJ = 1, N
            VTST = ABS( VR( JJ, J ) )
            IF( VTST.GT.VMX ) VMX = VTST             IF( AIMAG( VR( JJ, J ) ).EQ.ZERO .AND. ABS( REAL( VR( JJ, J ) ) ).GT.VRMX ) VRMX = ABS( REAL( VR( JJ, J ) ) )
   20    CONTINUE
         IF( VRMX / VMX.LT.ONE-TWO*ULP ) RESULT( 3 ) = ULPINV
   30 CONTINUE

      // Do Test (4)

      DO 50 J = 1, N
         TNRM = SCNRM2( N, VL( 1, J ), 1 )
         RESULT( 4 ) = MAX( RESULT( 4 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
         VMX = ZERO
         VRMX = ZERO
         DO 40 JJ = 1, N
            VTST = ABS( VL( JJ, J ) )
            IF( VTST.GT.VMX ) VMX = VTST             IF( AIMAG( VL( JJ, J ) ).EQ.ZERO .AND. ABS( REAL( VL( JJ, J ) ) ).GT.VRMX ) VRMX = ABS( REAL( VL( JJ, J ) ) )
   40    CONTINUE
         IF( VRMX / VMX.LT.ONE-TWO*ULP ) RESULT( 4 ) = ULPINV
   50 CONTINUE

      // Test for all options of computing condition numbers

      DO 200 ISENS = 1, ISENSM

         SENSE = SENS( ISENS )

         // Compute eigenvalues only, and test them

         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( BALANC, 'N', 'N', SENSE, N, H, LDA, W1, CDUM, 1, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO )
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX2', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX2', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5)

         DO 60 J = 1, N
            IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
   60    CONTINUE

         // Do Test (8)

         if ( .NOT.NOBAL ) {
            DO 70 J = 1, N
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
   70       CONTINUE
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         }

         // Do Test (9)

         if ( ISENS.EQ.2 .AND. N.GT.1 ) {
            DO 80 J = 1, N
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
   80       CONTINUE
         }

         // Compute eigenvalues and right eigenvectors, and test them

         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( BALANC, 'N', 'V', SENSE, N, H, LDA, W1, CDUM, 1, LRE, LDLRE, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO )
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX3', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX3', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5) again

         DO 90 J = 1, N
            IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
   90    CONTINUE

         // Do Test (6)

         DO 110 J = 1, N
            DO 100 JJ = 1, N
               IF( VR( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 6 ) = ULPINV
  100       CONTINUE
  110    CONTINUE

         // Do Test (8) again

         if ( .NOT.NOBAL ) {
            DO 120 J = 1, N
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
  120       CONTINUE
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         }

         // Do Test (9) again

         if ( ISENS.EQ.2 .AND. N.GT.1 ) {
            DO 130 J = 1, N
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
  130       CONTINUE
         }

         // Compute eigenvalues and left eigenvectors, and test them

         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( BALANC, 'V', 'N', SENSE, N, H, LDA, W1, LRE, LDLRE, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO )
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX4', IINFO, N, JTYPE, BALANC, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX4', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 190
         }

         // Do Test (5) again

         DO 140 J = 1, N
            IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
  140    CONTINUE

         // Do Test (7)

         DO 160 J = 1, N
            DO 150 JJ = 1, N
               IF( VL( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 7 ) = ULPINV
  150       CONTINUE
  160    CONTINUE

         // Do Test (8) again

         if ( .NOT.NOBAL ) {
            DO 170 J = 1, N
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
  170       CONTINUE
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         }

         // Do Test (9) again

         if ( ISENS.EQ.2 .AND. N.GT.1 ) {
            DO 180 J = 1, N
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
  180       CONTINUE
         }

  190    CONTINUE

  200 CONTINUE

      // If COMP, compare condition numbers to precomputed ones

      if ( COMP ) {
         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( 'N', 'V', 'V', 'B', N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO )
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX5', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         DO 220 I = 1, N - 1
            KMIN = I
            if ( ISRT.EQ.0 ) {
               VRIMIN = REAL( W( I ) )
            } else {
               VRIMIN = AIMAG( W( I ) )
            }
            DO 210 J = I + 1, N
               if ( ISRT.EQ.0 ) {
                  VRICMP = REAL( W( J ) )
               } else {
                  VRICMP = AIMAG( W( J ) )
               }
               if ( VRICMP.LT.VRIMIN ) {
                  KMIN = J
                  VRIMIN = VRICMP
               }
  210       CONTINUE
            CTMP = W( KMIN )
            W( KMIN ) = W( I )
            W( I ) = CTMP
            VRIMIN = RCONDE( KMIN )
            RCONDE( KMIN ) = RCONDE( I )
            RCONDE( I ) = VRIMIN
            VRIMIN = RCONDV( KMIN )
            RCONDV( KMIN ) = RCONDV( I )
            RCONDV( I ) = VRIMIN
  220    CONTINUE

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         RESULT( 10 ) = ZERO
         EPS = MAX( EPSIN, ULP )
         V = MAX( REAL( N )*EPS*ABNRM, SMLNUM )
         IF( ABNRM.EQ.ZERO ) V = ONE
         DO 230 I = 1, N
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
  230    CONTINUE

         // Compare condition numbers for eigenvalues
         // taking their condition numbers into account

         RESULT( 11 ) = ZERO
         DO 240 I = 1, N
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
  240    CONTINUE
  250    CONTINUE

      }

 9999 FORMAT( ' CGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' CGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', BALANC = ', A, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )

      RETURN

      // End of CGET23

      }
