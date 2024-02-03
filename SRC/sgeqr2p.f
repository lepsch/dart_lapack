      SUBROUTINE SGEQR2P( M, N, A, LDA, TAU, WORK, INFO );

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
      // EXTERNAL SLARF, SLARFGP, XERBLA
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
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SGEQR2P', -INFO );
         return;
      }

      K = MIN( M, N );

      for (I = 1; I <= K; I++) { // 10

         // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

         slarfgp(M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, TAU( I ) );
         if ( I < N ) {

            // Apply H(i) to A(i:m,i+1:n) from the left

            AII = A( I, I );
            A( I, I ) = ONE;
            slarf('Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK );
            A( I, I ) = AII;
         }
      } // 10
      return;

      // End of SGEQR2P

      }
