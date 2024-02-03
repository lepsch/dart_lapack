      SUBROUTINE SORG2L( M, N, K, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, II, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SSCAL, XERBLA
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
         xerbla('SORG2L', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      // Initialise columns 1:n-k to columns of the unit matrix

      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE

      DO 40 I = 1, K
         II = N - K + I

         // Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

         A( M-N+II, II ) = ONE
         slarf('Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, LDA, WORK );
         sscal(M-N+II-1, -TAU( I ), A( 1, II ), 1 );
         A( M-N+II, II ) = ONE - TAU( I )

         // Set A(m-k+i+1:m,n-k+i) to zero

         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN

      // End of SORG2L

      }
