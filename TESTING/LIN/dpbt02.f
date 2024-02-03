      SUBROUTINE DPBT02( UPLO, N, KD, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * );
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
      double             DASUM, DLAMCH, DLANSB;
      // EXTERNAL DASUM, DLAMCH, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSBMV
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

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSB( '1', UPLO, N, KD, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute  B - A*X

      for (J = 1; J <= NRHS; J++) { // 10
         dsbmv(UPLO, N, KD, -ONE, A, LDA, X( 1, J ), 1, ONE, B( 1, J ), 1 );
      } // 10

      // Compute the maximum over the number of right hand sides of
           // norm( B - A*X ) / ( norm(A) * norm(X) * EPS )

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 20
         BNORM = DASUM( N, B( 1, J ), 1 )
         XNORM = DASUM( N, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 20

      RETURN

      // End of DPBT02

      }
