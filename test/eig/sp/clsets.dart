      void clsets(final int M, final int P, final int N, final int A, final int AF, final int LDA, final int B, final int BF, final int LDB, final int C, final int CF, final int D, final int DF, final int X, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final int RESULT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, P, N;
      double               RESULT( 2 ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), C( * ), D( * ), CF( * ), DF( * ), WORK( LWORK ), X( * );

// ====================================================================

      int                INFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGGLSE, CLACPY, CGET02

      // Copy the matrices A and B to the arrays AF and BF,
      // and the vectors C and D to the arrays CF and DF,

      clacpy('Full', M, N, A, LDA, AF, LDA );
      clacpy('Full', P, N, B, LDB, BF, LDB );
      ccopy(M, C, 1, CF, 1 );
      ccopy(P, D, 1, DF, 1 );

      // Solve LSE problem

      cgglse(M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK, LWORK, INFO );

      // Test the residual for the solution of LSE

      // Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS

      ccopy(M, C, 1, CF, 1 );
      ccopy(P, D, 1, DF, 1 );
      cget02('No transpose', M, N, 1, A, LDA, X, N, CF, M, RWORK, RESULT( 1 ) );

      // Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS

      cget02('No transpose', P, N, 1, B, LDB, X, N, DF, P, RWORK, RESULT( 2 ) );

      }
