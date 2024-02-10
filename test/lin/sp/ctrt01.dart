      void ctrt01(UPLO, DIAG, N, final Matrix<double> A, final int LDA, final Matrix<double> AINV, final int LDAINV, RCOND, final Array<double> RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                LDA, LDAINV, N;
      double               RCOND, RESID;
      double               RWORK( * );
      Complex            A( LDA, * ), AINV( LDAINV, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANTR, SLAMCH;
      // EXTERNAL lsame, CLANTR, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANTR( '1', UPLO, DIAG, N, N, A, LDA, RWORK );
      AINVNM = CLANTR( '1', UPLO, DIAG, N, N, AINV, LDAINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Set the diagonal of AINV to 1 if AINV has unit diagonal.

      if ( lsame( DIAG, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            AINV[J][J] = ONE;
         } // 10
      }

      // Compute A * AINV, overwriting AINV.

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            ctrmv('Upper', 'No transpose', DIAG, J, A, LDA, AINV( 1, J ), 1 );
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 30
            ctrmv('Lower', 'No transpose', DIAG, N-J+1, A( J, J ), LDA, AINV( J, J ), 1 );
         } // 30
      }

      // Subtract 1 from each diagonal element to form A*AINV - I.

      for (J = 1; J <= N; J++) { // 40
         AINV[J][J] = AINV( J, J ) - ONE;
      } // 40

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, RWORK );

      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS;

      }
