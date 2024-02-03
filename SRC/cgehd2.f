      SUBROUTINE CGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, N;
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
      int                I;
      COMPLEX            ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      if ( N < 0 ) {
         INFO = -1
      } else if ( ILO < 1 || ILO > MAX( 1, N ) ) {
         INFO = -2
      } else if ( IHI < MIN( ILO, N ) || IHI > N ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      }
      if ( INFO != 0 ) {
         xerbla('CGEHD2', -INFO );
         RETURN
      }

      for (I = ILO; I <= IHI - 1; I++) { // 10

         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

         ALPHA = A( I+1, I )
         clarfg(IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) );
         A( I+1, I ) = ONE

         // Apply H(i) to A(1:ihi,i+1:ihi) from the right

         clarf('Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK );

         // Apply H(i)**H to A(i+1:ihi,i+1:n) from the left

         clarf('Left', IHI-I, N-I, A( I+1, I ), 1, CONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK );

         A( I+1, I ) = ALPHA
      } // 10

      RETURN

      // End of CGEHD2

      }
