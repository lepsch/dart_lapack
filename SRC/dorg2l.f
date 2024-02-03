      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, II, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || N > M ) {
         INFO = -2;
      } else if ( K < 0 || K > N ) {
         INFO = -3;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('DORG2L', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) RETURN;

      // Initialise columns 1:n-k to columns of the unit matrix

      for (J = 1; J <= N - K; J++) { // 20
         for (L = 1; L <= M; L++) { // 10
            A( L, J ) = ZERO;
         } // 10
         A( M-N+J, J ) = ONE;
      } // 20

      for (I = 1; I <= K; I++) { // 40
         II = N - K + I;

         // Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

         A( M-N+II, II ) = ONE;
         dlarf('Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, LDA, WORK );
         dscal(M-N+II-1, -TAU( I ), A( 1, II ), 1 );
         A( M-N+II, II ) = ONE - TAU( I );

         // Set A(m-k+i+1:m,n-k+i) to zero

         for (L = M - N + II + 1; L <= M; L++) { // 30
            A( L, II ) = ZERO;
         } // 30
      } // 40
      return;

      // End of DORG2L

      }
