      void sgeql2(M, N, A, LDA, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, K;
      REAL               AII;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SLARFG, XERBLA
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
         xerbla('SGEQL2', -INFO );
         return;
      }

      K = min( M, N );

      for (I = K; I >= 1; I--) { // 10

         // Generate elementary reflector H(i) to annihilate
         // A(1:m-k+i-1,n-k+i)

         slarfg(M-K+I, A( M-K+I, N-K+I ), A( 1, N-K+I ), 1, TAU( I ) );

         // Apply H(i) to A(1:m-k+i,1:n-k+i-1) from the left

         AII = A( M-K+I, N-K+I );
         A[M-K+I, N-K+I] = ONE;
         slarf('Left', M-K+I, N-K+I-1, A( 1, N-K+I ), 1, TAU( I ), A, LDA, WORK );
         A[M-K+I, N-K+I] = AII;
      } // 10
      return;
      }
