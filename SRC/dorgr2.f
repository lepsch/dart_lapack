      SUBROUTINE DORGR2( M, N, K, A, LDA, TAU, WORK, INFO )

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
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.M ) {
         INFO = -2
      } else if ( K.LT.0 .OR. K.GT.M ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('DORGR2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 ) RETURN

      if ( K.LT.M ) {

         // Initialise rows 1:m-k to rows of the unit matrix

         for (J = 1; J <= N; J++) { // 20
            DO 10 L = 1, M - K
               A( L, J ) = ZERO
            } // 10
            IF( J.GT.N-M .AND. J.LE.N-K ) A( M-N+J, J ) = ONE
         } // 20
      }

      for (I = 1; I <= K; I++) { // 40
         II = M - K + I

         // Apply H(i) to A(1:m-k+i,1:n-k+i) from the right

         A( II, N-M+II ) = ONE
         dlarf('Right', II-1, N-M+II, A( II, 1 ), LDA, TAU( I ), A, LDA, WORK );
         dscal(N-M+II-1, -TAU( I ), A( II, 1 ), LDA );
         A( II, N-M+II ) = ONE - TAU( I )

         // Set A(m-k+i,n-k+i+1:n) to zero

         DO 30 L = N - M + II + 1, N
            A( II, L ) = ZERO
         } // 30
      } // 40
      RETURN

      // End of DORGR2

      }
