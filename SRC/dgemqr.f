      SUBROUTINE DGEMQR( SIDE, TRANS, M, N, K, A, LDA, T, TSIZE, C, LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, LDA, M, N, K, TSIZE, LWORK, LDC
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), T( * ), C( LDC, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN, LQUERY
      INTEGER            MB, NB, LW, NBLCKS, MN, MINMNK, LWMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMQRT, DLAMTSQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      LQUERY  = ( LWORK.EQ.-1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'T' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
*
      MB = INT( T( 2 ) )
      NB = INT( T( 3 ) )
      IF( LEFT ) THEN
        LW = N * NB
        MN = M
      ELSE
        LW = MB * NB
        MN = N
      END IF
*
      MINMNK = MIN( M, N, K )
      IF( MINMNK.EQ.0 ) THEN
         LWMIN = 1
      ELSE
         LWMIN = MAX( 1, LW )
      END IF
*
      IF( ( MB.GT.K ) .AND. ( MN.GT.K ) ) THEN
        IF( MOD( MN - K, MB - K ).EQ.0 ) THEN
          NBLCKS = ( MN - K ) / ( MB - K )
        ELSE
          NBLCKS = ( MN - K ) / ( MB - K ) + 1
        END IF
      ELSE
        NBLCKS = 1
      END IF
*
      INFO = 0
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
        INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
        INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.MN ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, MN ) ) THEN
        INFO = -7
      ELSE IF( TSIZE.LT.5 ) THEN
        INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -11
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
        INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
        WORK( 1 ) = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGEMQR', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MINMNK.EQ.0 ) THEN
        RETURN
      END IF
*
      IF( ( LEFT .AND. M.LE.K ) .OR. ( RIGHT .AND. N.LE.K ) .OR. ( MB.LE.K ) .OR. ( MB.GE.MAX( M, N, K ) ) ) THEN         CALL DGEMQRT( SIDE, TRANS, M, N, K, NB, A, LDA, T( 6 ), NB, C, LDC, WORK, INFO )
      ELSE
        CALL DLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T( 6 ), NB, C, LDC, WORK, LWORK, INFO )
      END IF
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
*     End of DGEMQR
*
      END
