      double           FUNCTION ZQRT11( M, K, A, LDA, TAU, WORK, LWORK );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      int                INFO, J;
*     ..
*     .. External Functions ..
      double             DLAMCH, ZLANGE;
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASET, ZUNM2R
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX
*     ..
*     .. Local Arrays ..
      double             RDUMMY( 1 );
*     ..
*     .. Executable Statements ..
*
      ZQRT11 = ZERO
*
*     Test for sufficient workspace
*
      IF( LWORK.LT.M*M+M ) THEN
         CALL XERBLA( 'ZQRT11', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) RETURN
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, M )
*
*     Form Q
*
      CALL ZUNM2R( 'Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO )
*
*     Form Q'*Q
*
      CALL ZUNM2R( 'Left', 'Conjugate transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO )
*
      DO J = 1, M
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE
      END DO
*
      ZQRT11 = ZLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( DBLE( M )*DLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of ZQRT11
*
      END
