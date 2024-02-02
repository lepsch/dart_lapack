      SUBROUTINE SHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI,
     $                   VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL,
     $                   IFAILR, INFO )
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
      REAL               H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BOTHV, FROMQR, LEFTV, NOINIT, PAIR, RIGHTV
      INTEGER            I, IINFO, K, KL, KLN, KR, KSI, KSR, LDWORK
      REAL               BIGNUM, EPS3, HNORM, SMLNUM, ULP, UNFL, WKI,
     $                   WKR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      REAL               SLAMCH, SLANHS
      EXTERNAL           LSAME, SLAMCH, SLANHS, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAEIN, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
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
*     eigenvectors, and standardize the array SELECT.
*
      M = 0
      PAIR = .FALSE.
      DO 10 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
            SELECT( K ) = .FALSE.
         ELSE
            IF( WI( K ).EQ.ZERO ) THEN
               IF( SELECT( K ) )
     $            M = M + 1
            ELSE
               PAIR = .TRUE.
               IF( SELECT( K ) .OR. SELECT( K+1 ) ) THEN
                  SELECT( K ) = .TRUE.
                  M = M + 2
               END IF
            END IF
         END IF
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
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      ELSE IF( MM.LT.M ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SHSEIN', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Set machine-dependent constants.
*
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM
*
      LDWORK = N + 1
*
      KL = 1
      KLN = 0
      IF( FROMQR ) THEN
         KR = 0
      ELSE
         KR = N
      END IF
      KSR = 1
*
      DO 120 K = 1, N
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
                  IF( H( I, I-1 ).EQ.ZERO )
     $               GO TO 30
   20          CONTINUE
   30          CONTINUE
               KL = I
               IF( K.GT.KR ) THEN
                  DO 40 I = K, N - 1
                     IF( H( I+1, I ).EQ.ZERO )
     $                  GO TO 50
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
               HNORM = SLANHS( 'I', KR-KL+1, H( KL, KL ), LDH, WORK )
               IF( SISNAN( HNORM ) ) THEN
                  INFO = -6
                  RETURN
               ELSE IF( HNORM.GT.ZERO ) THEN
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
            WKR = WR( K )
            WKI = WI( K )
   60       CONTINUE
            DO 70 I = K - 1, KL, -1
               IF( SELECT( I ) .AND. ABS( WR( I )-WKR )+
     $             ABS( WI( I )-WKI ).LT.EPS3 ) THEN
                  WKR = WKR + EPS3
                  GO TO 60
               END IF
   70       CONTINUE
            WR( K ) = WKR
*
            PAIR = WKI.NE.ZERO
            IF( PAIR ) THEN
               KSI = KSR + 1
            ELSE
               KSI = KSR
            END IF
            IF( LEFTV ) THEN
*
*              Compute left eigenvector.
*
               CALL SLAEIN( .FALSE., NOINIT, N-KL+1, H( KL, KL ), LDH,
     $                      WKR, WKI, VL( KL, KSR ), VL( KL, KSI ),
     $                      WORK, LDWORK, WORK( N*N+N+1 ), EPS3, SMLNUM,
     $                      BIGNUM, IINFO )
               IF( IINFO.GT.0 ) THEN
                  IF( PAIR ) THEN
                     INFO = INFO + 2
                  ELSE
                     INFO = INFO + 1
                  END IF
                  IFAILL( KSR ) = K
                  IFAILL( KSI ) = K
               ELSE
                  IFAILL( KSR ) = 0
                  IFAILL( KSI ) = 0
               END IF
               DO 80 I = 1, KL - 1
                  VL( I, KSR ) = ZERO
   80          CONTINUE
               IF( PAIR ) THEN
                  DO 90 I = 1, KL - 1
                     VL( I, KSI ) = ZERO
   90             CONTINUE
               END IF
            END IF
            IF( RIGHTV ) THEN
*
*              Compute right eigenvector.
*
               CALL SLAEIN( .TRUE., NOINIT, KR, H, LDH, WKR, WKI,
     $                      VR( 1, KSR ), VR( 1, KSI ), WORK, LDWORK,
     $                      WORK( N*N+N+1 ), EPS3, SMLNUM, BIGNUM,
     $                      IINFO )
               IF( IINFO.GT.0 ) THEN
                  IF( PAIR ) THEN
                     INFO = INFO + 2
                  ELSE
                     INFO = INFO + 1
                  END IF
                  IFAILR( KSR ) = K
                  IFAILR( KSI ) = K
               ELSE
                  IFAILR( KSR ) = 0
                  IFAILR( KSI ) = 0
               END IF
               DO 100 I = KR + 1, N
                  VR( I, KSR ) = ZERO
  100          CONTINUE
               IF( PAIR ) THEN
                  DO 110 I = KR + 1, N
                     VR( I, KSI ) = ZERO
  110             CONTINUE
               END IF
            END IF
*
            IF( PAIR ) THEN
               KSR = KSR + 2
            ELSE
               KSR = KSR + 1
            END IF
         END IF
  120 CONTINUE
*
      RETURN
*
*     End of SHSEIN
*
      END