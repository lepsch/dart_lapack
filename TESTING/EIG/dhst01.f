      SUBROUTINE DHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK, LWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, LDA, LDH, LDQ, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), H( LDH, * ), Q( LDQ, * ), RESULT( 2 ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                LDWORK;
      double             ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DORT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RETURN
      }

      UNFL = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      OVFL = ONE / UNFL
      SMLNUM = UNFL*N / EPS

      // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

      // Copy A to WORK

      LDWORK = MAX( 1, N )
      dlacpy(' ', N, N, A, LDA, WORK, LDWORK );

      // Compute Q*H

      dgemm('No transpose', 'No transpose', N, N, N, ONE, Q, LDQ, H, LDH, ZERO, WORK( LDWORK*N+1 ), LDWORK );

      // Compute A - Q*H*Q'

      dgemm('No transpose', 'Transpose', N, N, N, -ONE, WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, ONE, WORK, LDWORK );

      ANORM = MAX( DLANGE( '1', N, N, A, LDA, WORK( LDWORK*N+1 ) ), UNFL )
      WNORM = DLANGE( '1', N, N, WORK, LDWORK, WORK( LDWORK*N+1 ) )

      // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

      RESULT( 1 ) = MIN( WNORM, ANORM ) / MAX( SMLNUM, ANORM*EPS ) / N

      // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

      dort01('Columns', N, N, Q, LDQ, WORK, LWORK, RESULT( 2 ) );

      RETURN

      // End of DHST01

      }
