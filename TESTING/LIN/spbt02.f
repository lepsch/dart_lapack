      SUBROUTINE SPBT02( UPLO, N, KD, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * )
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
      REAL               SASUM, SLAMCH, SLANSB
      // EXTERNAL SASUM, SLAMCH, SLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSB( '1', UPLO, N, KD, A, LDA, RWORK )
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute  B - A*X

      for (J = 1; J <= NRHS; J++) { // 10
         ssbmv(UPLO, N, KD, -ONE, A, LDA, X( 1, J ), 1, ONE, B( 1, J ), 1 );
      } // 10

      // Compute the maximum over the number of right hand sides of
           // norm( B - A*X ) / ( norm(A) * norm(X) * EPS )

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 20
         BNORM = SASUM( N, B( 1, J ), 1 )
         XNORM = SASUM( N, X( 1, J ), 1 )
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 20

      RETURN

      // End of SPBT02

      }
