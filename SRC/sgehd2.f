      SUBROUTINE SGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               AII
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEHD2', -INFO )
         RETURN
      END IF
*
      DO 10 I = ILO, IHI - 1
*
         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         CALL SLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, TAU( I ) )
         AII = A( I+1, I )
         A( I+1, I ) = ONE
*
         // Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL SLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK )
*
         // Apply H(i) to A(i+1:ihi,i+1:n) from the left
*
         CALL SLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ), A( I+1, I+1 ), LDA, WORK )
*
         A( I+1, I ) = AII
   10 CONTINUE
*
      RETURN
*
      // End of SGEHD2
*
      END
