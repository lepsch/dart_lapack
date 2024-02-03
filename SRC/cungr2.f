      SUBROUTINE CUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO
      const              ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, II, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACGV, CLARF, CSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX
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
         xerbla('CUNGR2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 ) RETURN

      if ( K.LT.M ) {

         // Initialise rows 1:m-k to rows of the unit matrix

         DO 20 J = 1, N
            DO 10 L = 1, M - K
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.N-M .AND. J.LE.N-K ) A( M-N+J, J ) = ONE
   20    CONTINUE
      }

      DO 40 I = 1, K
         II = M - K + I

         // Apply H(i)**H to A(1:m-k+i,1:n-k+i) from the right

         clacgv(N-M+II-1, A( II, 1 ), LDA );
         A( II, N-M+II ) = ONE
         clarf('Right', II-1, N-M+II, A( II, 1 ), LDA, CONJG( TAU( I ) ), A, LDA, WORK );
         cscal(N-M+II-1, -TAU( I ), A( II, 1 ), LDA );
         clacgv(N-M+II-1, A( II, 1 ), LDA );
         A( II, N-M+II ) = ONE - CONJG( TAU( I ) )

         // Set A(m-k+i,n-k+i+1:n) to zero

         DO 30 L = N - M + II + 1, N
            A( II, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN

      // End of CUNGR2

      }
