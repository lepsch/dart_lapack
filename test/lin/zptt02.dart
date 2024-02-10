      void zptt02(final int UPLO, final int N, final int NRHS, final int D, final int E, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      double             D( * );
      Complex         B( LDB, * ), E( * ), X( LDX, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DZASUM, ZLANHT;
      // EXTERNAL DLAMCH, DZASUM, ZLANHT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLAPTM

      // Quick return if possible

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Compute the 1-norm of the tridiagonal matrix A.

      ANORM = ZLANHT( '1', N, D, E );

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute B - A*X.

      zlaptm(UPLO, N, NRHS, -ONE, D, E, X, LDX, ONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
      //    norm(B - A*X) / ( norm(A) * norm(X) * EPS ).

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

      }
