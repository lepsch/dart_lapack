      SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX*16         ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARF, ZLARFG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) {
         INFO = -2
      } else if ( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGEHD2', -INFO );
         RETURN
      }

      for (I = ILO; I <= IHI - 1; I++) { // 10

         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

         ALPHA = A( I+1, I )
         zlarfg(IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) );
         A( I+1, I ) = ONE

         // Apply H(i) to A(1:ihi,i+1:ihi) from the right

         zlarf('Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK );

         // Apply H(i)**H to A(i+1:ihi,i+1:n) from the left

         zlarf('Left', IHI-I, N-I, A( I+1, I ), 1, DCONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK );

         A( I+1, I ) = ALPHA
      } // 10

      RETURN

      // End of ZGEHD2

      }
