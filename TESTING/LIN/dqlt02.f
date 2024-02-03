      SUBROUTINE DQLT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      double             ROGUE;
      PARAMETER          ( ROGUE = -1.0D+10 )
*     ..
*     .. Local Scalars ..
      int                INFO;
      double             ANORM, EPS, RESID;
*     ..
*     .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASET, DORGQL, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Scalars in Common ..
      String             SRNAMT;
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the last k columns of the factorization to the array Q
*
      CALL DLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      IF( K.LT.M ) CALL DLACPY( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, N-K+1 ), LDA )       IF( K.GT.1 ) CALL DLACPY( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, N-K+2 ), LDA )
*
*     Generate the last n columns of the matrix Q
*
      SRNAMT = 'DORGQL'
      CALL DORGQL( M, N, K, Q, LDA, TAU( N-K+1 ), WORK, LWORK, INFO )
*
*     Copy L(m-n+1:m,n-k+1:n)
*
      CALL DLASET( 'Full', N, K, ZERO, ZERO, L( M-N+1, N-K+1 ), LDA )
      CALL DLACPY( 'Lower', K, K, AF( M-K+1, N-K+1 ), LDA, L( M-K+1, N-K+1 ), LDA )
*
*     Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)
*
      CALL DGEMM( 'Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A( 1, N-K+1 ), LDA, ONE, L( M-N+1, N-K+1 ), LDA )
*
*     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', M, K, A( 1, N-K+1 ), LDA, RWORK )
      RESID = DLANGE( '1', N, K, L( M-N+1, N-K+1 ), LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, L, LDA )
      CALL DSYRK( 'Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, L, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', N, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of DQLT02
*
      END
