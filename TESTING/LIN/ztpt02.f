      SUBROUTINE ZTPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDX, B, LDB, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
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
      double             DLAMCH, DZASUM, ZLANTP;
      // EXTERNAL LSAME, DLAMCH, DZASUM, ZLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute the 1-norm of op(A).

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = ZLANTP( '1', UPLO, DIAG, N, AP, RWORK )
      } else {
         ANORM = ZLANTP( 'I', UPLO, DIAG, N, AP, RWORK )
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
         zcopy(N, X( 1, J ), 1, WORK, 1 );
         ztpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         zaxpy(N, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         BNORM = DZASUM( N, WORK, 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of ZTPT02

      }
