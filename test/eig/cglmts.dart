      void cglmts(N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N;
      double               RESULT;
      // ..
      // .. Array Arguments ..
      double               RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), D( * ), DF( * ), U( * ), WORK( LWORK ), X( * );

// ====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      Complex            CONE;
      const              CONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double               ANORM, BNORM, EPS, XNORM, YNORM, DNORM, UNFL;
      // ..
      // .. External Functions ..
      //- REAL               SCASUM, SLAMCH, CLANGE;
      // EXTERNAL SCASUM, SLAMCH, CLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACPY

      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );
      ANORM = max( CLANGE( '1', N, M, A, LDA, RWORK ), UNFL );
      BNORM = max( CLANGE( '1', N, P, B, LDB, RWORK ), UNFL );

      // Copy the matrices A and B to the arrays AF and BF,
      // and the vector D the array DF.

      clacpy('Full', N, M, A, LDA, AF, LDA );
      clacpy('Full', N, P, B, LDB, BF, LDB );
      ccopy(N, D, 1, DF, 1 );

      // Solve GLM problem

      cggglm(N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO );

      // Test the residual for the solution of LSE

                        // norm( d - A*x - B*u )
        // RESULT = -----------------------------------------
                 // (norm(A)+norm(B))*(norm(x)+norm(u))*EPS

      ccopy(N, D, 1, DF, 1 );
      cgemv('No transpose', N, M, -CONE, A, LDA, X, 1, CONE, DF, 1 );

      cgemv('No transpose', N, P, -CONE, B, LDB, U, 1, CONE, DF, 1 );

      DNORM = SCASUM( N, DF, 1 );
      XNORM = SCASUM( M, X, 1 ) + SCASUM( P, U, 1 );
      YNORM = ANORM + BNORM;

      if ( XNORM <= ZERO ) {
         RESULT = ZERO;
      } else {
         RESULT =  ( ( DNORM / YNORM ) / XNORM ) /EPS;
      }

      return;
      }