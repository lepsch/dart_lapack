      void zgeql2(M, N, final Matrix<double> A, final int LDA, TAU, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      Complex         A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      int                I, K;
      Complex         ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARF, ZLARFG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN

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
         xerbla('ZGEQL2', -INFO );
         return;
      }

      K = min( M, N );

      for (I = K; I >= 1; I--) { // 10

         // Generate elementary reflector H(i) to annihilate
         // A(1:m-k+i-1,n-k+i)

         ALPHA = A( M-K+I, N-K+I );
         zlarfg(M-K+I, ALPHA, A( 1, N-K+I ), 1, TAU( I ) );

         // Apply H(i)**H to A(1:m-k+i,1:n-k+i-1) from the left

         A[M-K+I][N-K+I] = ONE;
         zlarf('Left', M-K+I, N-K+I-1, A( 1, N-K+I ), 1, DCONJG( TAU( I ) ), A, LDA, WORK );
         A[M-K+I][N-K+I] = ALPHA;
      } // 10
      }
