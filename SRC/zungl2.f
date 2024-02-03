      SUBROUTINE ZUNGL2( M, N, K, A, LDA, TAU, WORK, INFO )

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
      int                I, J, L;
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
         xerbla('ZUNGL2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 ) RETURN

      if ( K.LT.M ) {

         // Initialise rows k+1:m to rows of the unit matrix

         for (J = 1; J <= N; J++) { // 20
            DO 10 L = K + 1, M
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.K .AND. J.LE.M ) A( J, J ) = ONE
   20    CONTINUE
      }

      DO 40 I = K, 1, -1

         // Apply H(i)**H to A(i:m,i:n) from the right

         if ( I.LT.N ) {
            zlacgv(N-I, A( I, I+1 ), LDA );
            if ( I.LT.M ) {
               A( I, I ) = ONE
               zlarf('Right', M-I, N-I+1, A( I, I ), LDA, DCONJG( TAU( I ) ), A( I+1, I ), LDA, WORK );
            }
            zscal(N-I, -TAU( I ), A( I, I+1 ), LDA );
            zlacgv(N-I, A( I, I+1 ), LDA );
         }
         A( I, I ) = ONE - DCONJG( TAU( I ) )

         // Set A(i,1:i-1) to zero

         DO 30 L = 1, I - 1
            A( I, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN

      // End of ZUNGL2

      }
