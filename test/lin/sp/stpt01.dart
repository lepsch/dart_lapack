      void stpt01(final int UPLO, final int DIAG, final int N, final int AP, final int AINVP, final int RCOND, final Array<double> _WORK_, final int RESID,) {
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                N;
      double               RCOND, RESID;
      double               AINVP( * ), AP( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               UNITD;
      int                J, JC;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANTP;
      // EXTERNAL lsame, SLAMCH, SLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPMV
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
      ANORM = SLANTP( '1', UPLO, DIAG, N, AP, WORK );
      AINVNM = SLANTP( '1', UPLO, DIAG, N, AINVP, WORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Compute A * AINV, overwriting AINV.

      UNITD = lsame( DIAG, 'U' );
      if ( lsame( UPLO, 'U' ) ) {
         JC = 1;
         for (J = 1; J <= N; J++) { // 10
            if (UNITD) AINVP( JC+J-1 ) = ONE;

            // Form the j-th column of A*AINV

            stpmv('Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 );

            // Subtract 1 from the diagonal

            AINVP[JC+J-1] = AINVP( JC+J-1 ) - ONE;
            JC = JC + J;
         } // 10
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 20
            if (UNITD) AINVP( JC ) = ONE;

            // Form the j-th column of A*AINV

            stpmv('Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 );

            // Subtract 1 from the diagonal

            AINVP[JC] = AINVP( JC ) - ONE;
            JC = JC + N - J + 1;
         } // 20
      }

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = SLANTP( '1', UPLO, 'Non-unit', N, AINVP, WORK );

      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS;

      }