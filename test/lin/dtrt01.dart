      void dtrt01(UPLO, DIAG, N, final Matrix<double> A, final int LDA, final Matrix<double> AINV, final int LDAINV, RCOND, final Array<double> _WORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                LDA, LDAINV, N;
      double             RCOND, RESID;
      double             A( LDA, * ), AINV( LDAINV, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANTR;
      // EXTERNAL lsame, DLAMCH, DLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      // Quick exit if N = 0

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = DLANTR( '1', UPLO, DIAG, N, N, A, LDA, WORK );
      AINVNM = DLANTR( '1', UPLO, DIAG, N, N, AINV, LDAINV, WORK );
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
            dtrmv('Upper', 'No transpose', DIAG, J, A, LDA, AINV( 1, J ), 1 );
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 30
            dtrmv('Lower', 'No transpose', DIAG, N-J+1, A( J, J ), LDA, AINV( J, J ), 1 );
         } // 30
      }

      // Subtract 1 from each diagonal element to form A*AINV - I.

      for (J = 1; J <= N; J++) { // 40
         AINV[J][J] = AINV( J, J ) - ONE;
      } // 40

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = DLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, WORK );

      RESID = ( ( RESID*RCOND ) / N.toDouble() ) / EPS;

      }
