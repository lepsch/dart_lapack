      SUBROUTINE DRQT01( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ROGUE
      PARAMETER          ( ROGUE = -1.0D+10 )
*     ..
*     .. Local Scalars ..
      int                INFO, MINMN
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGERQF, DLACPY, DLASET, DORGRQ, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      MINMN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'DGERQF'
      CALL DGERQF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy details of Q
*
      CALL DLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N ) CALL DLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M.GT.1 ) CALL DLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 ) CALL DLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA )
      END IF
*
*     Generate the n-by-n matrix Q
*
      SRNAMT = 'DORGRQ'
      CALL DORGRQ( N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL DLASET( 'Full', M, N, ZERO, ZERO, R, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 ) CALL DLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA )
      ELSE
         IF( M.GT.N .AND. N.GT.0 ) CALL DLACPY( 'Full', M-N, N, AF, LDA, R, LDA )          IF( N.GT.0 ) CALL DLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA )
      END IF
*
*     Compute R - A*Q'
*
      CALL DGEMM( 'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA )
*
*     Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
      RESID = DLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of DRQT01
*
      END
