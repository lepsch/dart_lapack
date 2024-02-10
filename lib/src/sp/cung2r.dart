      void cung2r(M, N, K, A, LDA, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, M, N;
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      int                I, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || N > M ) {
         INFO = -2;
      } else if ( K < 0 || K > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('CUNG2R', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) return;

      // Initialise columns k+1:n to columns of the unit matrix

      for (J = K + 1; J <= N; J++) { // 20
         for (L = 1; L <= M; L++) { // 10
            A[L][J] = ZERO;
         } // 10
         A[J][J] = ONE;
      } // 20

      for (I = K; I >= 1; I--) { // 40

         // Apply H(i) to A(i:m,i:n) from the left

         if ( I < N ) {
            A[I][I] = ONE;
            clarf('Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK );
         }
         if (I < M) cscal( M-I, -TAU( I ), A( I+1, I ), 1 );
         A[I][I] = ONE - TAU( I );

         // Set A(1:i-1,i) to zero

         for (L = 1; L <= I - 1; L++) { // 30
            A[L][I] = ZERO;
         } // 30
      } // 40
      }
