      void dpot06(UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID;
      double             A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE, NEGONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      const              NEGONE = -1.0 ;
      int                IFAIL, J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- int                idamax;
      //- double             DLAMCH, DLANSY;
      // EXTERNAL idamax, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS == 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = DLANSY( 'I', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute  B - A*X  and store in B.
      IFAIL=0;

      dsymm('Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = ABS(B(idamax( N, B( 1, J ), 1 ),J));
         XNORM = ABS(X(idamax( N, X( 1, J ), 1 ),J));
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      return;
      }
