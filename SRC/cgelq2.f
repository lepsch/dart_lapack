      SUBROUTINE CGELQ2( M, N, A, LDA, TAU, WORK, INFO )

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
      // EXTERNAL CLACGV, CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
         CALL XERBLA( 'CGELQ2', -INFO )
         RETURN
      }

      K = MIN( M, N )

      DO 10 I = 1, K

         // Generate elementary reflector H(i) to annihilate A(i,i+1:n)

         CALL CLACGV( N-I+1, A( I, I ), LDA )
         ALPHA = A( I, I )
         CALL CLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA, TAU( I ) )
         if ( I.LT.M ) {

            // Apply H(i) to A(i+1:m,i:n) from the right

            A( I, I ) = ONE
            CALL CLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), A( I+1, I ), LDA, WORK )
         }
         A( I, I ) = ALPHA
         CALL CLACGV( N-I+1, A( I, I ), LDA )
   10 CONTINUE
      RETURN

      // End of CGELQ2

      }
