      SUBROUTINE STPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDX, B, LDB, WORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SASUM, SLAMCH, SLANTP;
      // EXTERNAL LSAME, SASUM, SLAMCH, SLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, STPMV
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
         // norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

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

      return;

      // End of STPT02

      }
