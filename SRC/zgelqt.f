      SUBROUTINE ZGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDT, M, N, MB
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, IINFO, K
*     ..
*     .. External Subroutines ..
      EXTERNAL   ZGELQT3, ZLARFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( MB.LT.1 .OR. (MB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 ))THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.MB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGELQT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN
*
*     Blocked loop of length K
*
      DO I = 1, K,  MB
         IB = MIN( K-I+1, MB )
*
*     Compute the LQ factorization of the current block A(I:M,I:I+IB-1)
*
         CALL ZGELQT3( IB, N-I+1, A(I,I), LDA, T(1,I), LDT, IINFO )
         IF( I+IB.LE.M ) THEN
*
*     Update by applying H**T to A(I:M,I+IB:N) from the right
*
         CALL ZLARFB( 'R', 'N', 'F', 'R', M-I-IB+1, N-I+1, IB,
     $                   A( I, I ), LDA, T( 1, I ), LDT,
     $                   A( I+IB, I ), LDA, WORK , M-I-IB+1 )
         END IF
      END DO
      RETURN
*
*     End of ZGELQT
*
      END