      double zqrt11(final int M, final int K, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK,) {
  final A = A_.having();
  final WORK = WORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M;
      Complex         A( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
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
      double             RDUMMY( 1 );

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

      ZQRT11 = ZLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( M.toDouble()*dlamch( 'Epsilon' ) );

      }
