      void zsyt02(UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DZASUM, ZLANSY;
      // EXTERNAL DLAMCH, DZASUM, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute  B - A*X  (or  B - A'*X ) and store in B .

      zsymm('Left', UPLO, N, NRHS, -CONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = DZASUM( N, B( 1, J ), 1 );
         XNORM = DZASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      return;
      }