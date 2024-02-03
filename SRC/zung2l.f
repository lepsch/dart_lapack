      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )

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
      // EXTERNAL XERBLA, ZLARF, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 .OR. N.GT.M ) {
         INFO = -2
      } else if ( K.LT.0 .OR. K.GT.N ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('ZUNG2L', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      // Initialise columns 1:n-k to columns of the unit matrix

      DO 20 J = 1, N - K
         for (L = 1; L <= M; L++) { // 10
            A( L, J ) = ZERO
         } // 10
         A( M-N+J, J ) = ONE
      } // 20

      for (I = 1; I <= K; I++) { // 40
         II = N - K + I

         // Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

         A( M-N+II, II ) = ONE
         zlarf('Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, LDA, WORK );
         zscal(M-N+II-1, -TAU( I ), A( 1, II ), 1 );
         A( M-N+II, II ) = ONE - TAU( I )

         // Set A(m-k+i+1:m,n-k+i) to zero

         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
         } // 30
      } // 40
      RETURN

      // End of ZUNG2L

      }
