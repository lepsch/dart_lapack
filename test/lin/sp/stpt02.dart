      void stpt02(UPLO, TRANS, DIAG, N, NRHS, AP, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      double               RESID;
      double               AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SASUM, SLAMCH, SLANTP;
      // EXTERNAL lsame, SASUM, SLAMCH, SLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, STPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Compute the 1-norm of op(A).

      if ( lsame( TRANS, 'N' ) ) {
         ANORM = SLANTP( '1', UPLO, DIAG, N, AP, WORK );
      } else {
         ANORM = SLANTP( 'I', UPLO, DIAG, N, AP, WORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         scopy(N, X( 1, J ), 1, WORK, 1 );
         stpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         saxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );
         BNORM = SASUM( N, WORK, 1 );
         XNORM = SASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      }
