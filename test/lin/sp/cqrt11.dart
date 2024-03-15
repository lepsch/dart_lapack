      double cqrt11(final int M, final int K, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK,) {
  final A = A_.dim();
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M;
      Complex            A( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                INFO, J;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL
      double               RDUMMY( 1 );

      CQRT11 = ZERO;

      // Test for sufficient workspace

      if ( LWORK < M*M+M ) {
         xerbla('CQRT11', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0) return;

      claset('Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), WORK, M );

      // Form Q

      cunm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      // Form Q'*Q

      cunm2r('Left', 'Conjugate transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      for (J = 1; J <= M; J++) {
         WORK[( J-1 )*M+J] = WORK( ( J-1 )*M+J ) - ONE;
      }

      CQRT11 = CLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( REAL( M )*SLAMCH( 'Epsilon' ) );

      }