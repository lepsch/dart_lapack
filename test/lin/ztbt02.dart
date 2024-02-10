      void ztbt02(UPLO, TRANS, DIAG, N, KD, NRHS, final Matrix<double> AB, final int LDAB, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<double> _WORK, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      double             RESID;
      double             RWORK( * );
      Complex         AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DZASUM, ZLANTB;
      // EXTERNAL lsame, DLAMCH, DZASUM, ZLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZTBMV
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
         ANORM = ZLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, RWORK );
      } else {
         ANORM = ZLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, RWORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         zcopy(N, X( 1, J ), 1, WORK, 1 );
         ztbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
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
