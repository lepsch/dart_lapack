      SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             AII;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, DLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
         CALL XERBLA( 'DGEHD2', -INFO )
         RETURN
      }

      DO 10 I = ILO, IHI - 1

         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

         CALL DLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, TAU( I ) )
         AII = A( I+1, I )
         A( I+1, I ) = ONE

         // Apply H(i) to A(1:ihi,i+1:ihi) from the right

         CALL DLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK )

         // Apply H(i) to A(i+1:ihi,i+1:n) from the left

         CALL DLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ), A( I+1, I+1 ), LDA, WORK )

         A( I+1, I ) = AII
   10 CONTINUE

      RETURN

      // End of DGEHD2

      }
