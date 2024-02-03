      REAL             FUNCTION CQRT11( M, K, A, LDA, TAU, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL
      // ..
      // .. Local Arrays ..
      REAL               RDUMMY( 1 );
      // ..
      // .. Executable Statements ..

      CQRT11 = ZERO;

      // Test for sufficient workspace

      if ( LWORK < M*M+M ) {
         xerbla('CQRT11', 7 );
         RETURN;
      }

      // Quick return if possible

      if (M <= 0) RETURN;

      claset('Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), WORK, M );

      // Form Q

      cunm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      // Form Q'*Q

      cunm2r('Left', 'Conjugate transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      for (J = 1; J <= M; J++) {
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE;
      }

      CQRT11 = CLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( REAL( M )*SLAMCH( 'Epsilon' ) );

      RETURN;

      // End of CQRT11

      }
