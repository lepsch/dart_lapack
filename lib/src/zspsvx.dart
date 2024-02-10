      void zspsvx(FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
      int                IPIV( * );
      double             BERR( * ), FERR( * ), RWORK( * );
      Complex         AFP( * ), AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      bool               NOFACT;
      double             ANORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANSP;
      // EXTERNAL lsame, DLAMCH, ZLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZSPCON, ZSPRFS, ZSPTRF, ZSPTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      if ( !NOFACT && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('ZSPSVX', -INFO );
         return;
      }

      if ( NOFACT ) {

         // Compute the factorization A = U*D*U**T or A = L*D*L**T.

         zcopy(N*( N+1 ) / 2, AP, 1, AFP, 1 );
         zsptrf(UPLO, N, AFP, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = ZLANSP( 'I', UPLO, N, AP, RWORK );

      // Compute the reciprocal of the condition number of A.

      zspcon(UPLO, N, AFP, IPIV, ANORM, RCOND, WORK, INFO );

      // Compute the solution vectors X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zsptrs(UPLO, N, NRHS, AFP, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      zsprfs(UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < dlamch( 'Epsilon' ) ) INFO = N + 1;

      }
