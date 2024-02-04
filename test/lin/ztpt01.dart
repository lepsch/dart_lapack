      void ztpt01(UPLO, DIAG, N, AP, AINVP, RCOND, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         AINVP( * ), AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UNITD;
      int                J, JC;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANTP;
      // EXTERNAL lsame, DLAMCH, ZLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = ZLANTP( '1', UPLO, DIAG, N, AP, RWORK );
      AINVNM = ZLANTP( '1', UPLO, DIAG, N, AINVP, RWORK );
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

            // Form the j-th column of A*AINV.

            ztpmv('Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 );

            // Subtract 1 from the diagonal to form A*AINV - I.

            AINVP[JC+J-1] = AINVP( JC+J-1 ) - ONE;
            JC = JC + J;
         } // 10
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 20
            if (UNITD) AINVP( JC ) = ONE;

            // Form the j-th column of A*AINV.

            ztpmv('Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 );

            // Subtract 1 from the diagonal to form A*AINV - I.

            AINVP[JC] = AINVP( JC ) - ONE;
            JC = JC + N - J + 1;
         } // 20
      }

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = ZLANTP( '1', UPLO, 'Non-unit', N, AINVP, RWORK );

      RESID = ( ( RESID*RCOND ) / N.toDouble() ) / EPS;

      return;
      }