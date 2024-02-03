      SUBROUTINE ZUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      const              ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, II, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLARF, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < M ) {
         INFO = -2
      } else if ( K < 0 || K.GT.M ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      }
      if ( INFO != 0 ) {
         xerbla('ZUNGR2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (M.LE.0) RETURN;

      if ( K < M ) {

         // Initialise rows 1:m-k to rows of the unit matrix

         for (J = 1; J <= N; J++) { // 20
            for (L = 1; L <= M - K; L++) { // 10
               A( L, J ) = ZERO
            } // 10
            if (J.GT.N-M && J.LE.N-K) A( M-N+J, J ) = ONE;
         } // 20
      }

      for (I = 1; I <= K; I++) { // 40
         II = M - K + I

         // Apply H(i)**H to A(1:m-k+i,1:n-k+i) from the right

         zlacgv(N-M+II-1, A( II, 1 ), LDA );
         A( II, N-M+II ) = ONE
         zlarf('Right', II-1, N-M+II, A( II, 1 ), LDA, DCONJG( TAU( I ) ), A, LDA, WORK );
         zscal(N-M+II-1, -TAU( I ), A( II, 1 ), LDA );
         zlacgv(N-M+II-1, A( II, 1 ), LDA );
         A( II, N-M+II ) = ONE - DCONJG( TAU( I ) )

         // Set A(m-k+i,n-k+i+1:n) to zero

         for (L = N - M + II + 1; L <= N; L++) { // 30
            A( II, L ) = ZERO
         } // 30
      } // 40
      RETURN

      // End of ZUNGR2

      }
