      SUBROUTINE CGEQR2( M, N, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, K;
      COMPLEX            ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('CGEQR2', -INFO );
         RETURN
      }

      K = MIN( M, N )

      for (I = 1; I <= K; I++) { // 10

         // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

         clarfg(M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, TAU( I ) );
         if ( I.LT.N ) {

            // Apply H(i)**H to A(i:m,i+1:n) from the left

            ALPHA = A( I, I )
            A( I, I ) = ONE
            clarf('Left', M-I+1, N-I, A( I, I ), 1, CONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK );
            A( I, I ) = ALPHA
         }
      } // 10
      RETURN

      // End of CGEQR2

      }
