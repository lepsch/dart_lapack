      SUBROUTINE DPOT06( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, NEGONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      const              NEGONE = -1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                IFAIL, J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DLANSY;
      // EXTERNAL IDAMAX, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N.LE.0 .OR. NRHS.EQ.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSY( 'I', UPLO, N, A, LDA, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute  B - A*X  and store in B.
      IFAIL=0

      dsymm('Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = ABS(B(IDAMAX( N, B( 1, J ), 1 ),J))
         XNORM = ABS(X(IDAMAX( N, X( 1, J ), 1 ),J))
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of DPOT06

      }
