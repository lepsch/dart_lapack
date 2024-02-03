      void cgerq2(M, N, A, LDA, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, K;
      COMPLEX            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACGV, CLARF, CLARFG, XERBLA
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
         xerbla('CGERQ2', -INFO );
         return;
      }

      K = min( M, N );

      DO 10 I = K, 1, -1;

         // Generate elementary reflector H(i) to annihilate
         // A(m-k+i,1:n-k+i-1)

         clacgv(N-K+I, A( M-K+I, 1 ), LDA );
         ALPHA = A( M-K+I, N-K+I );
         clarfg(N-K+I, ALPHA, A( M-K+I, 1 ), LDA, TAU( I ) );

         // Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right

         A( M-K+I, N-K+I ) = ONE;
         clarf('Right', M-K+I-1, N-K+I, A( M-K+I, 1 ), LDA, TAU( I ), A, LDA, WORK );
         A( M-K+I, N-K+I ) = ALPHA;
         clacgv(N-K+I-1, A( M-K+I, 1 ), LDA );
      } // 10
      return;
      }
