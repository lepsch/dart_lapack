      void zgehd2(N, ILO, IHI, A, LDA, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      Complex         ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARF, ZLARFG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN
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
         xerbla('ZGEHD2', -INFO );
         return;
      }

      for (I = ILO; I <= IHI - 1; I++) { // 10

         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

         ALPHA = A( I+1, I );
         zlarfg(IHI-I, ALPHA, A( min( I+2, N ), I ), 1, TAU( I ) );
         A[I+1, I] = ONE;

         // Apply H(i) to A(1:ihi,i+1:ihi) from the right

         zlarf('Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK );

         // Apply H(i)**H to A(i+1:ihi,i+1:n) from the left

         zlarf('Left', IHI-I, N-I, A( I+1, I ), 1, DCONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK );

         A[I+1, I] = ALPHA;
      } // 10

      return;
      }
