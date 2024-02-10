      void ztrt02(UPLO, TRANS, DIAG, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, WORK, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID;
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DZASUM, ZLANTR;
      // EXTERNAL lsame, DLAMCH, DZASUM, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Compute the 1-norm of op(A).

      if ( lsame( TRANS, 'N' ) ) {
         ANORM = ZLANTR( '1', UPLO, DIAG, N, N, A, LDA, RWORK );
      } else {
         ANORM = ZLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS )

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         zcopy(N, X( 1, J ), 1, WORK, 1 );
         ztrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         zaxpy(N, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         BNORM = DZASUM( N, WORK, 1 );
         XNORM = DZASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      }
