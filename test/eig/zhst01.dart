      void zhst01(final int N, final int ILO, final int IHI, final Matrix<double> A_, final int LDA, final Matrix<double> H_, final int LDH, final Matrix<double> Q_, final int LDQ, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final A = A_.dim();
  final H = H_.dim();
  final Q = Q_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, LDA, LDH, LDQ, LWORK, N;
      double             RESULT( 2 ), RWORK( * );
      Complex         A( LDA, * ), H( LDH, * ), Q( LDQ, * ), WORK( LWORK );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                LDWORK;
      double             ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY, ZUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN

      // Quick return if possible

      if ( N <= 0 ) {
         RESULT[1] = ZERO;
         RESULT[2] = ZERO;
         return;
      }

      UNFL = dlamch( 'Safe minimum' );
      EPS = dlamch( 'Precision' );
      OVFL = ONE / UNFL;
      SMLNUM = UNFL*N / EPS;

      // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

      // Copy A to WORK

      LDWORK = max( 1, N );
      zlacpy(' ', N, N, A, LDA, WORK, LDWORK );

      // Compute Q*H

      zgemm('No transpose', 'No transpose', N, N, N, DCMPLX( ONE ), Q, LDQ, H, LDH, DCMPLX( ZERO ), WORK( LDWORK*N+1 ), LDWORK );

      // Compute A - Q*H*Q'

      zgemm('No transpose', 'Conjugate transpose', N, N, N, DCMPLX( -ONE ), WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, DCMPLX( ONE ), WORK, LDWORK );

      ANORM = max( ZLANGE( '1', N, N, A, LDA, RWORK ), UNFL );
      WNORM = ZLANGE( '1', N, N, WORK, LDWORK, RWORK );

      // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

      RESULT[1] = min( WNORM, ANORM ) / max( SMLNUM, ANORM*EPS ) / N;

      // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

      zunt01('Columns', N, N, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) );

      }
