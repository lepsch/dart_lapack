      SUBROUTINE ZPTT02( UPLO, N, NRHS, D, E, X, LDX, B, LDB, RESID );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             D( * );
      COMPLEX*16         B( LDB, * ), E( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DZASUM, ZLANHT;
      // EXTERNAL DLAMCH, DZASUM, ZLANHT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLAPTM
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         RESID = ZERO;
         RETURN;
      }

      // Compute the 1-norm of the tridiagonal matrix A.

      ANORM = ZLANHT( '1', N, D, E );

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         RETURN;
      }

      // Compute B - A*X.

      zlaptm(UPLO, N, NRHS, -ONE, D, E, X, LDX, ONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
         // norm(B - A*X) / ( norm(A) * norm(X) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = DZASUM( N, B( 1, J ), 1 );
         XNORM = DZASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      RETURN;

      // End of ZPTT02

      }
