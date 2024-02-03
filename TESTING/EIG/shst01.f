      SUBROUTINE SHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK, LWORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, LDA, LDH, LDQ, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), H( LDH, * ), Q( LDQ, * ), RESULT( 2 ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                LDWORK;
      REAL               ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SORT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         RESULT( 1 ) = ZERO;
         RESULT( 2 ) = ZERO;
         return;
      }

      UNFL = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );
      OVFL = ONE / UNFL;
      SMLNUM = UNFL*N / EPS;

      // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

      // Copy A to WORK

      LDWORK = max( 1, N );
      slacpy(' ', N, N, A, LDA, WORK, LDWORK );

      // Compute Q*H

      sgemm('No transpose', 'No transpose', N, N, N, ONE, Q, LDQ, H, LDH, ZERO, WORK( LDWORK*N+1 ), LDWORK );

      // Compute A - Q*H*Q'

      sgemm('No transpose', 'Transpose', N, N, N, -ONE, WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, ONE, WORK, LDWORK );

      ANORM = max( SLANGE( '1', N, N, A, LDA, WORK( LDWORK*N+1 ) ), UNFL );
      WNORM = SLANGE( '1', N, N, WORK, LDWORK, WORK( LDWORK*N+1 ) );

      // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

      RESULT( 1 ) = min( WNORM, ANORM ) / max( SMLNUM, ANORM*EPS ) / N;

      // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

      sort01('Columns', N, N, Q, LDQ, WORK, LWORK, RESULT( 2 ) );

      return;

      // End of SHST01

      }
