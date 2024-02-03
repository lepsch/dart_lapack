      SUBROUTINE CHPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               BERR( * ), FERR( * ), RWORK( * );
      COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOFACT;
      REAL               ANORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHP, SLAMCH;
      // EXTERNAL LSAME, CLANHP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHPCON, CHPRFS, CHPTRF, CHPTRS, CLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NOFACT = LSAME( FACT, 'N' );
      if ( !NOFACT && !LSAME( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CHPSVX', -INFO );
         RETURN;
      }

      if ( NOFACT ) {

         // Compute the factorization A = U*D*U**H or A = L*D*L**H.

         ccopy(N*( N+1 ) / 2, AP, 1, AFP, 1 );
         chptrf(UPLO, N, AFP, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            RETURN;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = CLANHP( 'I', UPLO, N, AP, RWORK );

      // Compute the reciprocal of the condition number of A.

      chpcon(UPLO, N, AFP, IPIV, ANORM, RCOND, WORK, INFO );

      // Compute the solution vectors X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      chptrs(UPLO, N, NRHS, AFP, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      chprfs(UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1;

      RETURN;

      // End of CHPSVX

      }
