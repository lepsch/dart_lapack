      void dtbt02(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDX, B, LDB, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DASUM, DLAMCH, DLANTB;
      // EXTERNAL LSAME, DASUM, DLAMCH, DLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DTBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Compute the 1-norm of op(A).

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = DLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, WORK );
      } else {
         ANORM = DLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, WORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         dcopy(N, X( 1, J ), 1, WORK, 1 );
         dtbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
         daxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );
         BNORM = DASUM( N, WORK, 1 );
         XNORM = DASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      return;

      // End of DTBT02

      }
