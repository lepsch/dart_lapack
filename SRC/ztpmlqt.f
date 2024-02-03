      SUBROUTINE ZTPMLQT( SIDE, TRANS, M, N, K, L, MB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT
*     ..
*     .. Array Arguments ..
      COMPLEX*16         V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, NB, LB, KF, LDAQ
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZTPRFB
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
      TRAN   = LSAME( TRANS, 'C' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF ( LEFT ) THEN
         LDAQ = MAX( 1, K )
      ELSE IF ( RIGHT ) THEN
         LDAQ = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( L.LT.0 .OR. L.GT.K ) THEN
         INFO = -6
      ELSE IF( MB.LT.1 .OR. (MB.GT.K .AND. K.GT.0) ) THEN
         INFO = -7
      ELSE IF( LDV.LT.K ) THEN
         INFO = -9
      ELSE IF( LDT.LT.MB ) THEN
         INFO = -11
      ELSE IF( LDA.LT.LDAQ ) THEN
         INFO = -13
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTPMLQT', -INFO )
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
            NB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = 0
            END IF
            CALL ZTPRFB( 'L', 'C', 'F', 'R', NB, N, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = NB-N+L-I+1
            END IF
            CALL ZTPRFB( 'R', 'N', 'F', 'R', M, NB, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      ELSE IF( LEFT .AND. TRAN ) THEN
*
         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = 0
            END IF
            CALL ZTPRFB( 'L', 'N', 'F', 'R', NB, N, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = NB-N+L-I+1
            END IF
            CALL ZTPRFB( 'R', 'C', 'F', 'R', M, NB, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      END IF
*
      RETURN
*
*     End of ZTPMLQT
*
      END
