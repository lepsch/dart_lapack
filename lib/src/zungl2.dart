      void zungl2(M, N, K, final Matrix<double> A, final int LDA, TAU, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, M, N;
      Complex         A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      int                I, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLARF, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( K < 0 || K > M ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('ZUNGL2', -INFO );
         return;
      }

      // Quick return if possible

      if (M <= 0) return;

      if ( K < M ) {

         // Initialise rows k+1:m to rows of the unit matrix

         for (J = 1; J <= N; J++) { // 20
            for (L = K + 1; L <= M; L++) { // 10
               A[L][J] = ZERO;
            } // 10
            if (J > K && J <= M) A( J, J ) = ONE;
         } // 20
      }

      for (I = K; I >= 1; I--) { // 40

         // Apply H(i)**H to A(i:m,i:n) from the right

         if ( I < N ) {
            zlacgv(N-I, A( I, I+1 ), LDA );
            if ( I < M ) {
               A[I][I] = ONE;
               zlarf('Right', M-I, N-I+1, A( I, I ), LDA, DCONJG( TAU( I ) ), A( I+1, I ), LDA, WORK );
            }
            zscal(N-I, -TAU( I ), A( I, I+1 ), LDA );
            zlacgv(N-I, A( I, I+1 ), LDA );
         }
         A[I][I] = ONE - DCONJG( TAU( I ) );

         // Set A(i,1:i-1) to zero

         for (L = 1; L <= I - 1; L++) { // 30
            A[I][L] = ZERO;
         } // 30
      } // 40
      }
