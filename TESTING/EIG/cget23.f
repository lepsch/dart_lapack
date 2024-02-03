      SUBROUTINE CGET23( COMP, ISRT, BALANC, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, H, W, W1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            COMP
      CHARACTER          BALANC
      int                INFO, ISRT, JTYPE, LDA, LDLRE, LDVL, LDVR, LWORK, N, NOUNIT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 )
      REAL               RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), RWORK( * ), SCALE( * ), SCALE1( * )
      COMPLEX            A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), W1( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
      REAL               EPSIN
      PARAMETER          ( EPSIN = 5.9605E-8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BALOK, NOBAL
      CHARACTER          SENSE
      int                I, IHI, IHI1, IINFO, ILO, ILO1, ISENS, ISENSM, J, JJ, KMIN       REAL               ABNRM, ABNRM1, EPS, SMLNUM, TNRM, TOL, TOLIN, ULP, ULPINV, V, VMAX, VMX, VRICMP, VRIMIN, VRMX, VTST
      COMPLEX            CTMP
*     ..
*     .. Local Arrays ..
      CHARACTER          SENS( 2 )
      REAL               RES( 2 )
      COMPLEX            CDUM( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SCNRM2, SLAMCH
      EXTERNAL           LSAME, SCNRM2, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEEVX, CGET22, CLACPY, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, MIN, REAL
*     ..
*     .. Data statements ..
      DATA               SENS / 'N', 'V' /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      NOBAL = LSAME( BALANC, 'N' )
      BALOK = NOBAL .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'B' )
      INFO = 0
      IF( ISRT.NE.0 .AND. ISRT.NE.1 ) THEN
         INFO = -2
      ELSE IF( .NOT.BALOK ) THEN
         INFO = -3
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -5
      ELSE IF( NOUNIT.LE.0 ) THEN
         INFO = -7
      ELSE IF( N.LT.0 ) THEN
         INFO = -8
      ELSE IF( LDA.LT.1 .OR. LDA.LT.N ) THEN
         INFO = -10
      ELSE IF( LDVL.LT.1 .OR. LDVL.LT.N ) THEN
         INFO = -15
      ELSE IF( LDVR.LT.1 .OR. LDVR.LT.N ) THEN
         INFO = -17
      ELSE IF( LDLRE.LT.1 .OR. LDLRE.LT.N ) THEN
         INFO = -19
      ELSE IF( LWORK.LT.2*N .OR. ( COMP .AND. LWORK.LT.2*N+N*N ) ) THEN
         INFO = -30
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGET23', -INFO )
         RETURN
      END IF
*
*     Quick return if nothing to do
*
      DO 10 I = 1, 11
         RESULT( I ) = -ONE
   10 CONTINUE
*
      IF( N.EQ.0 ) RETURN
*
*     More Important constants
*
      ULP = SLAMCH( 'Precision' )
      SMLNUM = SLAMCH( 'S' )
      ULPINV = ONE / ULP
*
*     Compute eigenvalues and eigenvectors, and test them
*
      IF( LWORK.GE.2*N+N*N ) THEN
         SENSE = 'B'
         ISENSM = 2
      ELSE
         SENSE = 'E'
         ISENSM = 1
      END IF
      CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
      CALL CGEEVX( BALANC, 'V', 'V', SENSE, N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO )
      IF( IINFO.NE.0 ) THEN
         RESULT( 1 ) = ULPINV
         IF( JTYPE.NE.22 ) THEN
            WRITE( NOUNIT, FMT = 9998 )'CGEEVX1', IINFO, N, JTYPE, BALANC, ISEED
         ELSE
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX1', IINFO, N, ISEED( 1 )
         END IF
         INFO = ABS( IINFO )
         RETURN
      END IF
*
*     Do Test (1)
*
      CALL CGET22( 'N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES )
      RESULT( 1 ) = RES( 1 )
*
*     Do Test (2)
*
      CALL CGET22( 'C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES )
      RESULT( 2 ) = RES( 1 )
*
*     Do Test (3)
*
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
*
*     Do Test (4)
*
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
*
*     Test for all options of computing condition numbers
*
      DO 200 ISENS = 1, ISENSM
*
         SENSE = SENS( ISENS )
*
*        Compute eigenvalues only, and test them
*
         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( BALANC, 'N', 'N', SENSE, N, H, LDA, W1, CDUM, 1, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX2', IINFO, N, JTYPE, BALANC, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX2', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 190
         END IF
*
*        Do Test (5)
*
         DO 60 J = 1, N
            IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
   60    CONTINUE
*
*        Do Test (8)
*
         IF( .NOT.NOBAL ) THEN
            DO 70 J = 1, N
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
   70       CONTINUE
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         END IF
*
*        Do Test (9)
*
         IF( ISENS.EQ.2 .AND. N.GT.1 ) THEN
            DO 80 J = 1, N
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
   80       CONTINUE
         END IF
*
*        Compute eigenvalues and right eigenvectors, and test them
*
         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( BALANC, 'N', 'V', SENSE, N, H, LDA, W1, CDUM, 1, LRE, LDLRE, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX3', IINFO, N, JTYPE, BALANC, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX3', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 190
         END IF
*
*        Do Test (5) again
*
         DO 90 J = 1, N
            IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
   90    CONTINUE
*
*        Do Test (6)
*
         DO 110 J = 1, N
            DO 100 JJ = 1, N
               IF( VR( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 6 ) = ULPINV
  100       CONTINUE
  110    CONTINUE
*
*        Do Test (8) again
*
         IF( .NOT.NOBAL ) THEN
            DO 120 J = 1, N
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
  120       CONTINUE
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         END IF
*
*        Do Test (9) again
*
         IF( ISENS.EQ.2 .AND. N.GT.1 ) THEN
            DO 130 J = 1, N
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
  130       CONTINUE
         END IF
*
*        Compute eigenvalues and left eigenvectors, and test them
*
         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( BALANC, 'V', 'N', SENSE, N, H, LDA, W1, LRE, LDLRE, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX4', IINFO, N, JTYPE, BALANC, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX4', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 190
         END IF
*
*        Do Test (5) again
*
         DO 140 J = 1, N
            IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
  140    CONTINUE
*
*        Do Test (7)
*
         DO 160 J = 1, N
            DO 150 JJ = 1, N
               IF( VL( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 7 ) = ULPINV
  150       CONTINUE
  160    CONTINUE
*
*        Do Test (8) again
*
         IF( .NOT.NOBAL ) THEN
            DO 170 J = 1, N
               IF( SCALE( J ).NE.SCALE1( J ) ) RESULT( 8 ) = ULPINV
  170       CONTINUE
            IF( ILO.NE.ILO1 ) RESULT( 8 ) = ULPINV             IF( IHI.NE.IHI1 ) RESULT( 8 ) = ULPINV             IF( ABNRM.NE.ABNRM1 ) RESULT( 8 ) = ULPINV
         END IF
*
*        Do Test (9) again
*
         IF( ISENS.EQ.2 .AND. N.GT.1 ) THEN
            DO 180 J = 1, N
               IF( RCONDV( J ).NE.RCNDV1( J ) ) RESULT( 9 ) = ULPINV
  180       CONTINUE
         END IF
*
  190    CONTINUE
*
  200 CONTINUE
*
*     If COMP, compare condition numbers to precomputed ones
*
      IF( COMP ) THEN
         CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL CGEEVX( 'N', 'V', 'V', 'B', N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX5', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
*        Sort eigenvalues and condition numbers lexicographically
*        to compare with inputs
*
         DO 220 I = 1, N - 1
            KMIN = I
            IF( ISRT.EQ.0 ) THEN
               VRIMIN = REAL( W( I ) )
            ELSE
               VRIMIN = AIMAG( W( I ) )
            END IF
            DO 210 J = I + 1, N
               IF( ISRT.EQ.0 ) THEN
                  VRICMP = REAL( W( J ) )
               ELSE
                  VRICMP = AIMAG( W( J ) )
               END IF
               IF( VRICMP.LT.VRIMIN ) THEN
                  KMIN = J
                  VRIMIN = VRICMP
               END IF
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
*
*        Compare condition numbers for eigenvectors
*        taking their condition numbers into account
*
         RESULT( 10 ) = ZERO
         EPS = MAX( EPSIN, ULP )
         V = MAX( REAL( N )*EPS*ABNRM, SMLNUM )
         IF( ABNRM.EQ.ZERO ) V = ONE
         DO 230 I = 1, N
            IF( V.GT.RCONDV( I )*RCONDE( I ) ) THEN
               TOL = RCONDV( I )
            ELSE
               TOL = V / RCONDE( I )
            END IF
            IF( V.GT.RCDVIN( I )*RCDEIN( I ) ) THEN
               TOLIN = RCDVIN( I )
            ELSE
               TOLIN = V / RCDEIN( I )
            END IF
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            IF( EPS*( RCDVIN( I )-TOLIN ).GT.RCONDV( I )+TOL ) THEN
               VMAX = ONE / EPS
            ELSE IF( RCDVIN( I )-TOLIN.GT.RCONDV( I )+TOL ) THEN
               VMAX = ( RCDVIN( I )-TOLIN ) / ( RCONDV( I )+TOL )
            ELSE IF( RCDVIN( I )+TOLIN.LT.EPS*( RCONDV( I )-TOL ) ) THEN
               VMAX = ONE / EPS
            ELSE IF( RCDVIN( I )+TOLIN.LT.RCONDV( I )-TOL ) THEN
               VMAX = ( RCONDV( I )-TOL ) / ( RCDVIN( I )+TOLIN )
            ELSE
               VMAX = ONE
            END IF
            RESULT( 10 ) = MAX( RESULT( 10 ), VMAX )
  230    CONTINUE
*
*        Compare condition numbers for eigenvalues
*        taking their condition numbers into account
*
         RESULT( 11 ) = ZERO
         DO 240 I = 1, N
            IF( V.GT.RCONDV( I ) ) THEN
               TOL = ONE
            ELSE
               TOL = V / RCONDV( I )
            END IF
            IF( V.GT.RCDVIN( I ) ) THEN
               TOLIN = ONE
            ELSE
               TOLIN = V / RCDVIN( I )
            END IF
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            IF( EPS*( RCDEIN( I )-TOLIN ).GT.RCONDE( I )+TOL ) THEN
               VMAX = ONE / EPS
            ELSE IF( RCDEIN( I )-TOLIN.GT.RCONDE( I )+TOL ) THEN
               VMAX = ( RCDEIN( I )-TOLIN ) / ( RCONDE( I )+TOL )
            ELSE IF( RCDEIN( I )+TOLIN.LT.EPS*( RCONDE( I )-TOL ) ) THEN
               VMAX = ONE / EPS
            ELSE IF( RCDEIN( I )+TOLIN.LT.RCONDE( I )-TOL ) THEN
               VMAX = ( RCONDE( I )-TOL ) / ( RCDEIN( I )+TOLIN )
            ELSE
               VMAX = ONE
            END IF
            RESULT( 11 ) = MAX( RESULT( 11 ), VMAX )
  240    CONTINUE
  250    CONTINUE
*
      END IF
*
 9999 FORMAT( ' CGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' CGET23: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', BALANC = ', A, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of CGET23
*
      END
