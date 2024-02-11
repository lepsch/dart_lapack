      void stbt02(final int UPLO, final int TRANS, final int DIAG, final int N, final int KD, final int NRHS, final Matrix<double> AB, final int LDAB, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<double> _WORK, final int RESID,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      double               RESID;
      double               AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SASUM, SLAMCH, SLANTB;
      // EXTERNAL lsame, SASUM, SLAMCH, SLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, STBMV
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
         ANORM = SLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, WORK );
      } else {
         ANORM = SLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, WORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         scopy(N, X( 1, J ), 1, WORK, 1 );
         stbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
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
