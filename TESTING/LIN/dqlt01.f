      SUBROUTINE DQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      double             ROGUE;
      PARAMETER          ( ROGUE = -1.0D+10 )
      // ..
      // .. Local Scalars ..
      int                INFO, MINMN;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DGEQLF, DLACPY, DLASET, DORGQL, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )

      // Copy the matrix A to the array AF.

      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )

      // Factorize the matrix A in the array AF.

      SRNAMT = 'DGEQLF'
      CALL DGEQLF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )

      // Copy details of Q

      CALL DLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      IF( M.GE.N ) THEN
         IF( N.LT.M .AND. N.GT.0 ) CALL DLACPY( 'Full', M-N, N, AF, LDA, Q( 1, M-N+1 ), LDA )          IF( N.GT.1 ) CALL DLACPY( 'Upper', N-1, N-1, AF( M-N+1, 2 ), LDA, Q( M-N+1, M-N+2 ), LDA )
      ELSE
         IF( M.GT.1 ) CALL DLACPY( 'Upper', M-1, M-1, AF( 1, N-M+2 ), LDA, Q( 1, 2 ), LDA )
      END IF

      // Generate the m-by-m matrix Q

      SRNAMT = 'DORGQL'
      CALL DORGQL( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy L

      CALL DLASET( 'Full', M, N, ZERO, ZERO, L, LDA )
      IF( M.GE.N ) THEN
         IF( N.GT.0 ) CALL DLACPY( 'Lower', N, N, AF( M-N+1, 1 ), LDA, L( M-N+1, 1 ), LDA )
      ELSE
         IF( N.GT.M .AND. M.GT.0 ) CALL DLACPY( 'Full', M, N-M, AF, LDA, L, LDA )          IF( M.GT.0 ) CALL DLACPY( 'Lower', M, M, AF( 1, N-M+1 ), LDA, L( 1, N-M+1 ), LDA )
      END IF

      // Compute L - Q'*A

      CALL DGEMM( 'Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A, LDA, ONE, L, LDA )

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
      RESID = DLANGE( '1', M, N, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q'*Q

      CALL DLASET( 'Full', M, M, ZERO, ONE, L, LDA )
      CALL DSYRK( 'Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, L, LDA )

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = DLANSY( '1', 'Upper', M, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS

      RETURN

      // End of DQLT01

      END
