      SUBROUTINE DGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, C, LDC, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDC, M, N, MB, LDT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, LDWORK, KF, Q
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DLARFB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'T' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF( LEFT ) THEN
         LDWORK = MAX( 1, N )
         Q = M
      ELSE IF ( RIGHT ) THEN
         LDWORK = MAX( 1, M )
         Q = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.Q ) THEN
         INFO = -5
      ELSE IF( MB.LT.1 .OR. (MB.GT.K .AND. K.GT.0)) THEN
         INFO = -6
      ELSE IF( LDV.LT.MAX( 1, K ) ) THEN
          INFO = -8
      ELSE IF( LDT.LT.MB ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEMLQT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. NOTRAN ) THEN
*
         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            CALL DLARFB( 'L', 'T', 'F', 'R', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            CALL DLARFB( 'R', 'N', 'F', 'R', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( LEFT .AND. TRAN ) THEN
*
         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            CALL DLARFB( 'L', 'N', 'F', 'R', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            CALL DLARFB( 'R', 'T', 'F', 'R', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      END IF
*
      RETURN
*
*     End of DGEMLQT
*
      END
