      SUBROUTINE ZGEQL2( M, N, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, K;
      COMPLEX*16         ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARF, ZLARFG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQL2', -INFO )
         RETURN
      END IF

      K = MIN( M, N )

      DO 10 I = K, 1, -1

         // Generate elementary reflector H(i) to annihilate
         // A(1:m-k+i-1,n-k+i)

         ALPHA = A( M-K+I, N-K+I )
         CALL ZLARFG( M-K+I, ALPHA, A( 1, N-K+I ), 1, TAU( I ) )

         // Apply H(i)**H to A(1:m-k+i,1:n-k+i-1) from the left

         A( M-K+I, N-K+I ) = ONE
         CALL ZLARF( 'Left', M-K+I, N-K+I-1, A( 1, N-K+I ), 1, DCONJG( TAU( I ) ), A, LDA, WORK )
         A( M-K+I, N-K+I ) = ALPHA
   10 CONTINUE
      RETURN

      // End of ZGEQL2

      END
