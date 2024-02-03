      void zgelq2(M, N, A, LDA, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, K;
      Complex         ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLARF, ZLARFG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZGELQ2', -INFO );
         return;
      }

      K = min( M, N );

      for (I = 1; I <= K; I++) { // 10

         // Generate elementary reflector H(i) to annihilate A(i,i+1:n)

         zlacgv(N-I+1, A( I, I ), LDA );
         ALPHA = A( I, I );
         zlarfg(N-I+1, ALPHA, A( I, min( I+1, N ) ), LDA, TAU( I ) );
         if ( I < M ) {

            // Apply H(i) to A(i+1:m,i:n) from the right

            A( I, I ) = ONE;
            zlarf('Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), A( I+1, I ), LDA, WORK );
         }
         A( I, I ) = ALPHA;
         zlacgv(N-I+1, A( I, I ), LDA );
      } // 10
      return;

      // End of ZGELQ2

      }
