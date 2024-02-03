      REAL             FUNCTION SQRT11( M, K, A, LDA, TAU, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( LWORK )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Local Arrays ..
      REAL               RDUMMY( 1 )
      // ..
      // .. Executable Statements ..
*
      SQRT11 = ZERO
*
      // Test for sufficient workspace
*
      IF( LWORK.LT.M*M+M ) THEN
         CALL XERBLA( 'SQRT11', 7 )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.LE.0 ) RETURN
*
      CALL SLASET( 'Full', M, M, ZERO, ONE, WORK, M )
*
      // Form Q
*
      CALL SORM2R( 'Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO )
*
      // Form Q'*Q
*
      CALL SORM2R( 'Left', 'Transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO )
*
      DO J = 1, M
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE
      END DO
*
      SQRT11 = SLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( REAL( M )*SLAMCH( 'Epsilon' ) )
*
      RETURN
*
      // End of SQRT11
*
      END
