      SUBROUTINE SGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N;
      REAL               RESULT
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), RWORK( * ), D( * ), DF( * ), U( * ), WORK( LWORK ), X( * )
*
*  ====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                INFO;
      REAL               ANORM, BNORM, EPS, XNORM, YNORM, DNORM, UNFL
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      // EXTERNAL SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      // EXTERNAL SLACPY
*
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      ANORM = MAX( SLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( SLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vector D the array DF.
*
      CALL SLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL SLACPY( 'Full', N, P, B, LDB, BF, LDB )
      CALL SCOPY( N, D, 1, DF, 1 )
*
*     Solve GLM problem
*
      CALL SGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO )
*
*     Test the residual for the solution of LSE
*
*                       norm( d - A*x - B*u )
*       RESULT = -----------------------------------------
*                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
      CALL SCOPY( N, D, 1, DF, 1 )
      CALL SGEMV( 'No transpose', N, M, -ONE, A, LDA, X, 1, ONE, DF, 1 )
*
      CALL SGEMV( 'No transpose', N, P, -ONE, B, LDB, U, 1, ONE, DF, 1 )
*
      DNORM = SASUM( N, DF, 1 )
      XNORM = SASUM( M, X, 1 ) + SASUM( P, U, 1 )
      YNORM = ANORM + BNORM
*
      IF( XNORM.LE.ZERO ) THEN
         RESULT = ZERO
      ELSE
         RESULT =  ( ( DNORM / YNORM ) / XNORM ) /EPS
      END IF
*
      RETURN
*
*     End of SGLMTS
*
      END
