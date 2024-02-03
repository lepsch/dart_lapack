      SUBROUTINE DLQT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             ROGUE;
      const              ROGUE = -1.0D+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASET, DORGLQ, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )

      // Copy the first k rows of the factorization to the array Q

      dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      dlacpy('Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the first n columns of the matrix Q

      SRNAMT = 'DORGLQ'
      dorglq(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L(1:k,1:m)

      dlaset('Full', K, M, ZERO, ZERO, L, LDA );
      dlacpy('Lower', K, M, AF, LDA, L, LDA );

      // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'

      dgemm('No transpose', 'Transpose', K, M, N, -ONE, A, LDA, Q, LDA, ONE, L, LDA );

      // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .

      ANORM = DLANGE( '1', K, N, A, LDA, RWORK )
      RESID = DLANGE( '1', K, M, L, LDA, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute I - Q*Q'

      dlaset('Full', M, M, ZERO, ONE, L, LDA );
      dsyrk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = DLANSY( '1', 'Upper', M, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of DLQT02

      }
