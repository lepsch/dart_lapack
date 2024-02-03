      void ctrt01(UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                LDA, LDAINV, N;
      REAL               RCOND, RESID;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            A( LDA, * ), AINV( LDAINV, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANTR, SLAMCH;
      // EXTERNAL LSAME, CLANTR, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

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

      if ( LSAME( DIAG, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            AINV( J, J ) = ONE;
         } // 10
      }

      // Compute A * AINV, overwriting AINV.

      if ( LSAME( UPLO, 'U' ) ) {
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
         AINV( J, J ) = AINV( J, J ) - ONE;
      } // 40

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, RWORK );

      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS;

      return;
      }
