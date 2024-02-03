      SUBROUTINE SLQT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      REAL               ROGUE
      const              ROGUE = -1.0E+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      REAL               ANORM, EPS, RESID
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASET, SORGLQ, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' )

      // Copy the first k rows of the factorization to the array Q

      CALL SLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL SLACPY( 'Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA )

      // Generate the first n columns of the matrix Q

      SRNAMT = 'SORGLQ'
      CALL SORGLQ( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy L(1:k,1:m)

      CALL SLASET( 'Full', K, M, ZERO, ZERO, L, LDA )
      CALL SLACPY( 'Lower', K, M, AF, LDA, L, LDA )

      // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'

      CALL SGEMM( 'No transpose', 'Transpose', K, M, N, -ONE, A, LDA, Q, LDA, ONE, L, LDA )

      // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .

      ANORM = SLANGE( '1', K, N, A, LDA, RWORK )
      RESID = SLANGE( '1', K, M, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q*Q'

      CALL SLASET( 'Full', M, M, ZERO, ONE, L, LDA )
      CALL SSYRK( 'Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA )

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = SLANSY( '1', 'Upper', M, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of SLQT02

      }
