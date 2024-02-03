      REAL             FUNCTION CQRT11( M, K, A, LDA, TAU, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      int                INFO, J;
*     ..
*     .. External Functions ..
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASET, CUNM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL
*     ..
*     .. Local Arrays ..
      REAL               RDUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      CQRT11 = ZERO
*
*     Test for sufficient workspace
*
      IF( LWORK.LT.M*M+M ) THEN
         CALL XERBLA( 'CQRT11', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) RETURN
*
      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), WORK, M )
*
*     Form Q
*
      CALL CUNM2R( 'Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO )
*
*     Form Q'*Q
*
      CALL CUNM2R( 'Left', 'Conjugate transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO )
*
      DO J = 1, M
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE
      END DO
*
      CQRT11 = CLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( REAL( M )*SLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of CQRT11
*
      END
