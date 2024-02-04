      double zqrt11(M, K, A, LDA, TAU, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLASET, ZUNM2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Local Arrays ..
      double             RDUMMY( 1 );
      // ..
      // .. Executable Statements ..

      ZQRT11 = ZERO;

      // Test for sufficient workspace

      if ( LWORK < M*M+M ) {
         xerbla('ZQRT11', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0) return;

      zlaset('Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, M );

      // Form Q

      zunm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      // Form Q'*Q

      zunm2r('Left', 'Conjugate transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      for (J = 1; J <= M; J++) {
         WORK[( J-1 )*M+J] = WORK( ( J-1 )*M+J ) - ONE;
      }

      ZQRT11 = ZLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( DBLE( M )*DLAMCH( 'Epsilon' ) );

      return;
      }
