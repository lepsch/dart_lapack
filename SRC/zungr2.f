      SUBROUTINE ZUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGR2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      IF( K.LT.M ) THEN
*
*        Initialise rows 1:m-k to rows of the unit matrix
*
         DO 20 J = 1, N
            DO 10 L = 1, M - K
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.N-M .AND. J.LE.N-K )
     $         A( M-N+J, J ) = ONE
   20    CONTINUE
      END IF
*
      DO 40 I = 1, K
         II = M - K + I
*
*        Apply H(i)**H to A(1:m-k+i,1:n-k+i) from the right
*
         CALL ZLACGV( N-M+II-1, A( II, 1 ), LDA )
         A( II, N-M+II ) = ONE
         CALL ZLARF( 'Right', II-1, N-M+II, A( II, 1 ), LDA,
     $               DCONJG( TAU( I ) ), A, LDA, WORK )
         CALL ZSCAL( N-M+II-1, -TAU( I ), A( II, 1 ), LDA )
         CALL ZLACGV( N-M+II-1, A( II, 1 ), LDA )
         A( II, N-M+II ) = ONE - DCONJG( TAU( I ) )
*
*        Set A(m-k+i,n-k+i+1:n) to zero
*
         DO 30 L = N - M + II + 1, N
            A( II, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNGR2
*
      END