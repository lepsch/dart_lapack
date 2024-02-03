      void zglmts(N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, P;
      double             RESULT;
      // ..
      // .. Array Arguments ..

// ====================================================================

      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), D( * ), DF( * ), U( * ), WORK( LWORK ), X( * );
      // ..
      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      COMPLEX*16         CONE;
      const              CONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..

      // EXTERNAL ZCOPY, ZGEMV, ZGGGLM, ZLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' );
      UNFL = DLAMCH( 'Safe minimum' );
      ANORM = max( ZLANGE( '1', N, M, A, LDA, RWORK ), UNFL );
      BNORM = max( ZLANGE( '1', N, P, B, LDB, RWORK ), UNFL );

      // Copy the matrices A and B to the arrays AF and BF,
      // and the vector D the array DF.

      zlacpy('Full', N, M, A, LDA, AF, LDA );
      zlacpy('Full', N, P, B, LDB, BF, LDB );
      zcopy(N, D, 1, DF, 1 );

      // Solve GLM problem

      zggglm(N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO );

      // Test the residual for the solution of LSE

                        // norm( d - A*x - B*u )
        // RESULT = -----------------------------------------
                 // (norm(A)+norm(B))*(norm(x)+norm(u))*EPS

      zcopy(N, D, 1, DF, 1 );
      zgemv('No transpose', N, M, -CONE, A, LDA, X, 1, CONE, DF, 1 );

      zgemv('No transpose', N, P, -CONE, B, LDB, U, 1, CONE, DF, 1 );

      DNORM = DZASUM( N, DF, 1 );
      XNORM = DZASUM( M, X, 1 ) + DZASUM( P, U, 1 );
      YNORM = ANORM + BNORM;

      if ( XNORM <= ZERO ) {
         RESULT = ZERO;
      } else {
         RESULT = ( ( DNORM / YNORM ) / XNORM ) / EPS;
      }

      return;

      // End of ZGLMTS

      }
