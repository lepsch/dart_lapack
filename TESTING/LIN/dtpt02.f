      SUBROUTINE DTPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDX, B, LDB, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DASUM, DLAMCH, DLANTP;
      // EXTERNAL LSAME, DASUM, DLAMCH, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N.LE.0 || NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute the 1-norm of op(A).

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = DLANTP( '1', UPLO, DIAG, N, AP, WORK )
      } else {
         ANORM = DLANTP( 'I', UPLO, DIAG, N, AP, WORK )
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 10
         dcopy(N, X( 1, J ), 1, WORK, 1 );
         dtpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         daxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );
         BNORM = DASUM( N, WORK, 1 )
         XNORM = DASUM( N, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of DTPT02

      }
