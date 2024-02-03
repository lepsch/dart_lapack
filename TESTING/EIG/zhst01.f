      SUBROUTINE ZHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK, LWORK, RWORK, RESULT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, LDA, LDH, LDQ, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( 2 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), H( LDH, * ), Q( LDQ, * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                LDWORK;
      double             ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY, ZUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         RESULT( 1 ) = ZERO;
         RESULT( 2 ) = ZERO;
         return;
      }

      UNFL = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Precision' );
      OVFL = ONE / UNFL;
      SMLNUM = UNFL*N / EPS;

      // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

      // Copy A to WORK

      LDWORK = MAX( 1, N );
      zlacpy(' ', N, N, A, LDA, WORK, LDWORK );

      // Compute Q*H

      zgemm('No transpose', 'No transpose', N, N, N, DCMPLX( ONE ), Q, LDQ, H, LDH, DCMPLX( ZERO ), WORK( LDWORK*N+1 ), LDWORK );

      // Compute A - Q*H*Q'

      zgemm('No transpose', 'Conjugate transpose', N, N, N, DCMPLX( -ONE ), WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, DCMPLX( ONE ), WORK, LDWORK );

      ANORM = MAX( ZLANGE( '1', N, N, A, LDA, RWORK ), UNFL );
      WNORM = ZLANGE( '1', N, N, WORK, LDWORK, RWORK );

      // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

      RESULT( 1 ) = MIN( WNORM, ANORM ) / MAX( SMLNUM, ANORM*EPS ) / N;

      // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

      zunt01('Columns', N, N, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) );

      return;

      // End of ZHST01

      }
