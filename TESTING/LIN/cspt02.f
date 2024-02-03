      SUBROUTINE CSPT02( UPLO, N, NRHS, A, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, BNORM, EPS, XNORM
      // ..
      // .. External Functions ..
      REAL               CLANSP, SCASUM, SLAMCH
      // EXTERNAL CLANSP, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANSP( '1', UPLO, N, A, RWORK )
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute  B - A*X  for the matrix of right hand sides B.

      for (J = 1; J <= NRHS; J++) { // 10
         cspmv(UPLO, N, -CONE, A, X( 1, J ), 1, CONE, B( 1, J ), 1 );
      } // 10

      // Compute the maximum over the number of right hand sides of
         // norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 20
         BNORM = SCASUM( N, B( 1, J ), 1 )
         XNORM = SCASUM( N, X( 1, J ), 1 )
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM/ANORM )/XNORM )/EPS )
         }
      } // 20

      RETURN

      // End of CSPT02

      }
