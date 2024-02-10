      void zlsets(M, P, N, A, AF, LDA, B, BF, LDB, C, CF, D, DF, X, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final int RESULT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, N, P;

// ====================================================================

      double             RESULT( 2 ), RWORK( * );
      Complex         A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), C( * ), CF( * ), D( * ), DF( * ), WORK( LWORK ), X( * );
      int                INFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZGET02, ZGGLSE, ZLACPY

      // Copy the matrices A and B to the arrays AF and BF,
      // and the vectors C and D to the arrays CF and DF,

      zlacpy('Full', M, N, A, LDA, AF, LDA );
      zlacpy('Full', P, N, B, LDB, BF, LDB );
      zcopy(M, C, 1, CF, 1 );
      zcopy(P, D, 1, DF, 1 );

      // Solve LSE problem

      zgglse(M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK, LWORK, INFO );

      // Test the residual for the solution of LSE

      // Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS

      zcopy(M, C, 1, CF, 1 );
      zcopy(P, D, 1, DF, 1 );
      zget02('No transpose', M, N, 1, A, LDA, X, N, CF, M, RWORK, RESULT( 1 ) );

      // Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS

      zget02('No transpose', P, N, 1, B, LDB, X, N, DF, P, RWORK, RESULT( 2 ) );

      }
