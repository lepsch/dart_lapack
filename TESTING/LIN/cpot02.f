      SUBROUTINE CPOT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, BNORM, EPS, XNORM
      // ..
      // .. External Functions ..
      REAL               CLANHE, SCASUM, SLAMCH
      // EXTERNAL CLANHE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute  B - A*X

      chemm('Left', UPLO, N, NRHS, -CONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = SCASUM( N, B( 1, J ), 1 )
         XNORM = SCASUM( N, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM/ANORM )/XNORM )/EPS )
         }
      } // 10

      RETURN

      // End of CPOT02

      }
