      void ctpt01(UPLO, DIAG, N, AP, AINVP, RCOND, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                N;
      REAL               RCOND, RESID;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            AINVP( * ), AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UNITD;
      int                J, JC;
      REAL               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- REAL               CLANTP, SLAMCH;
      // EXTERNAL LSAME, CLANTP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANTP( '1', UPLO, DIAG, N, AP, RWORK );
      AINVNM = CLANTP( '1', UPLO, DIAG, N, AINVP, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Compute A * AINV, overwriting AINV.

      UNITD = LSAME( DIAG, 'U' );
      if ( LSAME( UPLO, 'U' ) ) {
         JC = 1;
         for (J = 1; J <= N; J++) { // 10
            if (UNITD) AINVP( JC+J-1 ) = ONE;

            // Form the j-th column of A*AINV.

            ctpmv('Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 );

            // Subtract 1 from the diagonal to form A*AINV - I.

            AINVP( JC+J-1 ) = AINVP( JC+J-1 ) - ONE;
            JC = JC + J;
         } // 10
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 20
            if (UNITD) AINVP( JC ) = ONE;

            // Form the j-th column of A*AINV.

            ctpmv('Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 );

            // Subtract 1 from the diagonal to form A*AINV - I.

            AINVP( JC ) = AINVP( JC ) - ONE;
            JC = JC + N - J + 1;
         } // 20
      }

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANTP( '1', UPLO, 'Non-unit', N, AINVP, RWORK );

      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS;

      return;
      }
