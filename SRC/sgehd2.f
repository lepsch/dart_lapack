      void sgehd2(N, ILO, IHI, A, LDA, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               AII;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -2;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SGEHD2', -INFO );
         return;
      }

      for (I = ILO; I <= IHI - 1; I++) { // 10

         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

         slarfg(IHI-I, A( I+1, I ), A( min( I+2, N ), I ), 1, TAU( I ) );
         AII = A( I+1, I );
         A( I+1, I ) = ONE;

         // Apply H(i) to A(1:ihi,i+1:ihi) from the right

         slarf('Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK );

         // Apply H(i) to A(i+1:ihi,i+1:n) from the left

         slarf('Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ), A( I+1, I+1 ), LDA, WORK );

         A( I+1, I ) = AII;
      } // 10

      return;

      // End of SGEHD2

      }
