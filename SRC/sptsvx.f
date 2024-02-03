      SUBROUTINE SPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), E( * ), EF( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

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
      REAL               SLAMCH, SLANST;
      // EXTERNAL LSAME, SLAMCH, SLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLACPY, SPTCON, SPTRFS, SPTTRF, SPTTRS, XERBLA
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
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('SPTSVX', -INFO );
         return;
      }

      if ( NOFACT ) {

         // Compute the L*D*L**T (or U**T*D*U) factorization of A.

         scopy(N, D, 1, DF, 1 );
         if (N > 1) CALL SCOPY( N-1, E, 1, EF, 1 );
         spttrf(N, DF, EF, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = SLANST( '1', N, D, E );

      // Compute the reciprocal of the condition number of A.

      sptcon(N, DF, EF, ANORM, RCOND, WORK, INFO );

      // Compute the solution vectors X.

      slacpy('Full', N, NRHS, B, LDB, X, LDX );
      spttrs(N, NRHS, DF, EF, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      sptrfs(N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1;

      return;

      // End of SPTSVX

      }
