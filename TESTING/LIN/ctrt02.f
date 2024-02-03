      SUBROUTINE CTRT02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDX, B, LDB, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * )
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
      REAL               CLANTR, SCASUM, SLAMCH
      // EXTERNAL LSAME, CLANTR, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute the 1-norm of op(A).

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = CLANTR( '1', UPLO, DIAG, N, N, A, LDA, RWORK )
      } else {
         ANORM = CLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK )
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS )

      RESID = ZERO
      DO 10 J = 1, NRHS
         ccopy(N, X( 1, J ), 1, WORK, 1 );
         ctrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         caxpy(N, CMPLX( -ONE ), B( 1, J ), 1, WORK, 1 );
         BNORM = SCASUM( N, WORK, 1 )
         XNORM = SCASUM( N, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
   10 CONTINUE

      RETURN

      // End of CTRT02

      }
