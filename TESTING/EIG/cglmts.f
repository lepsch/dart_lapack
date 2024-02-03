      SUBROUTINE CGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, WORK, LWORK, RWORK, RESULT )
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
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), D( * ), DF( * ), U( * ), WORK( LWORK ), X( * )
*
*  ====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                INFO;
      REAL               ANORM, BNORM, EPS, XNORM, YNORM, DNORM, UNFL
*     ..
*     .. External Functions ..
      REAL               SCASUM, SLAMCH, CLANGE
      EXTERNAL           SCASUM, SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY
*
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      ANORM = MAX( CLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( CLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vector D the array DF.
*
      CALL CLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', N, P, B, LDB, BF, LDB )
      CALL CCOPY( N, D, 1, DF, 1 )
*
*     Solve GLM problem
*
      CALL CGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO )
*
*     Test the residual for the solution of LSE
*
*                       norm( d - A*x - B*u )
*       RESULT = -----------------------------------------
*                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
      CALL CCOPY( N, D, 1, DF, 1 )
      CALL CGEMV( 'No transpose', N, M, -CONE, A, LDA, X, 1, CONE, DF, 1 )
*
      CALL CGEMV( 'No transpose', N, P, -CONE, B, LDB, U, 1, CONE, DF, 1 )
*
      DNORM = SCASUM( N, DF, 1 )
      XNORM = SCASUM( M, X, 1 ) + SCASUM( P, U, 1 )
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
*     End of CGLMTS
*
      END
