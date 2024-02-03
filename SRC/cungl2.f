      SUBROUTINE CUNGL2( M, N, K, A, LDA, TAU, WORK, INFO )

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
      int                I, J, L;
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
         CALL XERBLA( 'CUNGL2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 ) RETURN

      if ( K.LT.M ) {

         // Initialise rows k+1:m to rows of the unit matrix

         DO 20 J = 1, N
            DO 10 L = K + 1, M
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.K .AND. J.LE.M ) A( J, J ) = ONE
   20    CONTINUE
      }

      DO 40 I = K, 1, -1

         // Apply H(i)**H to A(i:m,i:n) from the right

         if ( I.LT.N ) {
            CALL CLACGV( N-I, A( I, I+1 ), LDA )
            if ( I.LT.M ) {
               A( I, I ) = ONE
               CALL CLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, CONJG( TAU( I ) ), A( I+1, I ), LDA, WORK )
            }
            CALL CSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
            CALL CLACGV( N-I, A( I, I+1 ), LDA )
         }
         A( I, I ) = ONE - CONJG( TAU( I ) )

         // Set A(i,1:i-1,i) to zero

         DO 30 L = 1, I - 1
            A( I, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN

      // End of CUNGL2

      }
