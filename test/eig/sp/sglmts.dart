      void sglmts(N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, final Array<double> WORK, final int LWORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, P, N;
      double               RESULT;
      double               A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), RWORK( * ), D( * ), DF( * ), U( * ), WORK( LWORK ), X( * );

// ====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                INFO;
      double               ANORM, BNORM, EPS, XNORM, YNORM, DNORM, UNFL;
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SLAMCH, SLANGE;
      // EXTERNAL SASUM, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY

      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );
      ANORM = max( SLANGE( '1', N, M, A, LDA, RWORK ), UNFL );
      BNORM = max( SLANGE( '1', N, P, B, LDB, RWORK ), UNFL );

      // Copy the matrices A and B to the arrays AF and BF,
      // and the vector D the array DF.

      slacpy('Full', N, M, A, LDA, AF, LDA );
      slacpy('Full', N, P, B, LDB, BF, LDB );
      scopy(N, D, 1, DF, 1 );

      // Solve GLM problem

      sggglm(N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO );

      // Test the residual for the solution of LSE

                        // norm( d - A*x - B*u )
        // RESULT = -----------------------------------------
        //          (norm(A)+norm(B))*(norm(x)+norm(u))*EPS

      scopy(N, D, 1, DF, 1 );
      sgemv('No transpose', N, M, -ONE, A, LDA, X, 1, ONE, DF, 1 );

      sgemv('No transpose', N, P, -ONE, B, LDB, U, 1, ONE, DF, 1 );

      DNORM = SASUM( N, DF, 1 );
      XNORM = SASUM( M, X, 1 ) + SASUM( P, U, 1 );
      YNORM = ANORM + BNORM;

      if ( XNORM <= ZERO ) {
         RESULT = ZERO;
      } else {
         RESULT =  ( ( DNORM / YNORM ) / XNORM ) /EPS;
      }

      }
