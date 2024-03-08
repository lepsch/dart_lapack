      void dtpt02(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final int AP, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final Array<double> WORK_, final int RESID,) {
  final X = X_.having();
  final B = B_.having();
  final WORK = WORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      double             AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DASUM, DLAMCH, DLANTP;
      // EXTERNAL lsame, DASUM, DLAMCH, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DTPMV
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
         ANORM = DLANTP( '1', UPLO, DIAG, N, AP, WORK );
      } else {
         ANORM = DLANTP( 'I', UPLO, DIAG, N, AP, WORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         dcopy(N, X( 1, J ), 1, WORK, 1 );
         dtpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         daxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );
         BNORM = dasum( N, WORK, 1 );
         XNORM = dasum( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      }
