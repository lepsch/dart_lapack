      SUBROUTINE DGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, P;
      double             RESULT;
      // ..
      // .. Array Arguments ..

*  ====================================================================

      double             A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), D( * ), DF( * ), RWORK( * ), U( * ), WORK( LWORK ), X( * );
      // ..
      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..

      // EXTERNAL DCOPY, DGEMV, DGGGLM, DLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      ANORM = MAX( DLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( DLANGE( '1', N, P, B, LDB, RWORK ), UNFL )

      // Copy the matrices A and B to the arrays AF and BF,
      // and the vector D the array DF.

      dlacpy('Full', N, M, A, LDA, AF, LDA );
      dlacpy('Full', N, P, B, LDB, BF, LDB );
      dcopy(N, D, 1, DF, 1 );

      // Solve GLM problem

      dggglm(N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO );

      // Test the residual for the solution of LSE

                        // norm( d - A*x - B*u )
        // RESULT = -----------------------------------------
                 // (norm(A)+norm(B))*(norm(x)+norm(u))*EPS

      dcopy(N, D, 1, DF, 1 );
      dgemv('No transpose', N, M, -ONE, A, LDA, X, 1, ONE, DF, 1 );

      dgemv('No transpose', N, P, -ONE, B, LDB, U, 1, ONE, DF, 1 );

      DNORM = DASUM( N, DF, 1 )
      XNORM = DASUM( M, X, 1 ) + DASUM( P, U, 1 )
      YNORM = ANORM + BNORM

      if ( XNORM.LE.ZERO ) {
         RESULT = ZERO
      } else {
         RESULT = ( ( DNORM / YNORM ) / XNORM ) / EPS
      }

      RETURN

      // End of DGLMTS

      }
