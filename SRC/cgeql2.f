      SUBROUTINE CGEQL2( M, N, A, LDA, TAU, WORK, INFO );

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
      // EXTERNAL CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CGEQL2', -INFO );
         return;
      }

      K = MIN( M, N );

      DO 10 I = K, 1, -1;

         // Generate elementary reflector H(i) to annihilate
         // A(1:m-k+i-1,n-k+i)

         ALPHA = A( M-K+I, N-K+I );
         clarfg(M-K+I, ALPHA, A( 1, N-K+I ), 1, TAU( I ) );

         // Apply H(i)**H to A(1:m-k+i,1:n-k+i-1) from the left

         A( M-K+I, N-K+I ) = ONE;
         clarf('Left', M-K+I, N-K+I-1, A( 1, N-K+I ), 1, CONJG( TAU( I ) ), A, LDA, WORK );
         A( M-K+I, N-K+I ) = ALPHA;
      } // 10
      return;

      // End of CGEQL2

      }
