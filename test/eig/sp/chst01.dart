      void chst01(N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, LDA, LDH, LDQ, LWORK, N;
      double               RESULT( 2 ), RWORK( * );
      Complex            A( LDA, * ), H( LDH, * ), Q( LDQ, * ), WORK( LWORK );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                LDWORK;
      double               ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN

      // Quick return if possible

      if ( N <= 0 ) {
         RESULT[1] = ZERO;
         RESULT[2] = ZERO;
         return;
      }

      UNFL = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );
      OVFL = ONE / UNFL;
      SMLNUM = UNFL*N / EPS;

      // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

      // Copy A to WORK

      LDWORK = max( 1, N );
      clacpy(' ', N, N, A, LDA, WORK, LDWORK );

      // Compute Q*H

      cgemm('No transpose', 'No transpose', N, N, N, CMPLX( ONE ), Q, LDQ, H, LDH, CMPLX( ZERO ), WORK( LDWORK*N+1 ), LDWORK );

      // Compute A - Q*H*Q'

      cgemm('No transpose', 'Conjugate transpose', N, N, N, CMPLX( -ONE ), WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, CMPLX( ONE ), WORK, LDWORK );

      ANORM = max( CLANGE( '1', N, N, A, LDA, RWORK ), UNFL );
      WNORM = CLANGE( '1', N, N, WORK, LDWORK, RWORK );

      // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

      RESULT[1] = min( WNORM, ANORM ) / max( SMLNUM, ANORM*EPS ) / N;

      // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

      cunt01('Columns', N, N, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) );

      return;
      }
