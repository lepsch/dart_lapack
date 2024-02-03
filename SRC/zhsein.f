      SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, IFAILR, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EIGSRC, INITV, SIDE
      INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IFAILL( * ), IFAILR( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   RZERO
      PARAMETER          ( RZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BOTHV, FROMQR, LEFTV, NOINIT, RIGHTV
      INTEGER            I, IINFO, K, KL, KLN, KR, KS, LDWORK
      DOUBLE PRECISION   EPS3, HNORM, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM, WK
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           LSAME, DLAMCH, ZLANHS, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLAEIN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters.
*
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
*
      FROMQR = LSAME( EIGSRC, 'Q' )
*
      NOINIT = LSAME( INITV, 'N' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      M = 0
      DO 10 K = 1, N
         IF( SELECT( K ) ) M = M + 1
   10 CONTINUE
*
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.FROMQR .AND. .NOT.LSAME( EIGSRC, 'N' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOINIT .AND. .NOT.LSAME( INITV, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -10
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -12
      ELSE IF( MM.LT.M ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHSEIN', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) RETURN
*
*     Set machine-dependent constants.
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
*
      LDWORK = N
*
      KL = 1
      KLN = 0
      IF( FROMQR ) THEN
         KR = 0
      ELSE
         KR = N
      END IF
      KS = 1
*
      DO 100 K = 1, N
         IF( SELECT( K ) ) THEN
*
*           Compute eigenvector(s) corresponding to W(K).
*
            IF( FROMQR ) THEN
*
*              If affiliation of eigenvalues is known, check whether
*              the matrix splits.
*
*              Determine KL and KR such that 1 <= KL <= K <= KR <= N
*              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
*              KR = N).
*
*              Then inverse iteration can be performed with the
*              submatrix H(KL:N,KL:N) for a left eigenvector, and with
*              the submatrix H(1:KR,1:KR) for a right eigenvector.
*
               DO 20 I = K, KL + 1, -1
                  IF( H( I, I-1 ).EQ.ZERO ) GO TO 30
   20          CONTINUE
   30          CONTINUE
               KL = I
               IF( K.GT.KR ) THEN
                  DO 40 I = K, N - 1
                     IF( H( I+1, I ).EQ.ZERO ) GO TO 50
   40             CONTINUE
   50             CONTINUE
                  KR = I
               END IF
            END IF
*
            IF( KL.NE.KLN ) THEN
               KLN = KL
*
*              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
*              has not ben computed before.
*
               HNORM = ZLANHS( 'I', KR-KL+1, H( KL, KL ), LDH, RWORK )
               IF( DISNAN( HNORM ) ) THEN
                  INFO = -6
                  RETURN
               ELSE IF( HNORM.GT.RZERO ) THEN
                  EPS3 = HNORM*ULP
               ELSE
                  EPS3 = SMLNUM
               END IF
            END IF
*
*           Perturb eigenvalue if it is close to any previous
*           selected eigenvalues affiliated to the submatrix
*           H(KL:KR,KL:KR). Close roots are modified by EPS3.
*
            WK = W( K )
   60       CONTINUE
            DO 70 I = K - 1, KL, -1
               IF( SELECT( I ) .AND. CABS1( W( I )-WK ).LT.EPS3 ) THEN
                  WK = WK + EPS3
                  GO TO 60
               END IF
   70       CONTINUE
            W( K ) = WK
*
            IF( LEFTV ) THEN
*
*              Compute left eigenvector.
*
               CALL ZLAEIN( .FALSE., NOINIT, N-KL+1, H( KL, KL ), LDH, WK, VL( KL, KS ), WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO )
               IF( IINFO.GT.0 ) THEN
                  INFO = INFO + 1
                  IFAILL( KS ) = K
               ELSE
                  IFAILL( KS ) = 0
               END IF
               DO 80 I = 1, KL - 1
                  VL( I, KS ) = ZERO
   80          CONTINUE
            END IF
            IF( RIGHTV ) THEN
*
*              Compute right eigenvector.
*
               CALL ZLAEIN( .TRUE., NOINIT, KR, H, LDH, WK, VR( 1, KS ), WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO )
               IF( IINFO.GT.0 ) THEN
                  INFO = INFO + 1
                  IFAILR( KS ) = K
               ELSE
                  IFAILR( KS ) = 0
               END IF
               DO 90 I = KR + 1, N
                  VR( I, KS ) = ZERO
   90          CONTINUE
            END IF
            KS = KS + 1
         END IF
  100 CONTINUE
*
      RETURN
*
*     End of ZHSEIN
*
      END
