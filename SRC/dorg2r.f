      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )

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
      int                I, J, L;
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
      } else if ( N.LT.0 .OR. N.GT.M ) {
         INFO = -2
      } else if ( K.LT.0 .OR. K.GT.N ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('DORG2R', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      // Initialise columns k+1:n to columns of the unit matrix

      DO 20 J = K + 1, N
         for (L = 1; L <= M; L++) { // 10
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE

      DO 40 I = K, 1, -1

         // Apply H(i) to A(i:m,i:n) from the left

         if ( I.LT.N ) {
            A( I, I ) = ONE
            dlarf('Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK );
         }
         IF( I.LT.M ) CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )

         // Set A(1:i-1,i) to zero

         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN

      // End of DORG2R

      }
