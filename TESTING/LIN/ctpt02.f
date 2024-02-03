      SUBROUTINE CTPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDX, B, LDB, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, BNORM, EPS, XNORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANTP, SCASUM, SLAMCH
      // EXTERNAL LSAME, CLANTP, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute the 1-norm of op(A).

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = CLANTP( '1', UPLO, DIAG, N, AP, RWORK )
      } else {
         ANORM = CLANTP( 'I', UPLO, DIAG, N, AP, RWORK )
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 10
         ccopy(N, X( 1, J ), 1, WORK, 1 );
         ctpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         caxpy(N, CMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         BNORM = SCASUM( N, WORK, 1 )
         XNORM = SCASUM( N, X( 1, J ), 1 )
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of CTPT02

      }
