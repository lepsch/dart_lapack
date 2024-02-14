      void spot03(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AINV_, final int LDAINV, final Matrix<double> WORK_, final int LDWORK, final Array<double> RWORK_, final int RCOND, final int RESID,) {
  final A = A_.dim();
  final AINV = AINV_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAINV, LDWORK, N;
      double               RCOND, RESID;
      double               A( LDA, * ), AINV( LDAINV, * ), RWORK( * ), WORK( LDWORK, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL lsame, SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK );
      AINVNM = SLANSY( '1', UPLO, N, AINV, LDAINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Expand AINV into a full matrix and call SSYMM to multiply
      // AINV on the left by A.

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J - 1; I++) { // 10
               AINV[J][I] = AINV( I, J );
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J + 1; I <= N; I++) { // 30
               AINV[J][I] = AINV( I, J );
            } // 30
         } // 40
      }
      ssymm('Left', UPLO, N, N, -ONE, A, LDA, AINV, LDAINV, ZERO, WORK, LDWORK );

      // Add the identity matrix to WORK .

      for (I = 1; I <= N; I++) { // 50
         WORK[I][I] = WORK( I, I ) + ONE;
      } // 50

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = SLANGE( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND ) / EPS ) / REAL( N );

      }
