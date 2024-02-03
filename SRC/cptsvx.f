      SUBROUTINE CPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      REAL               BERR( * ), D( * ), DF( * ), FERR( * ), RWORK( * )       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOFACT;
      REAL               ANORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHT, SLAMCH
      // EXTERNAL LSAME, CLANHT, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CPTCON, CPTRFS, CPTTRF, CPTTRS, SCOPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      if ( !NOFACT && !LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -11
      }
      if ( INFO != 0 ) {
         xerbla('CPTSVX', -INFO );
         RETURN
      }

      if ( NOFACT ) {

         // Compute the L*D*L**H (or U**H*D*U) factorization of A.

         scopy(N, D, 1, DF, 1 );
         if (N > 1) CALL CCOPY( N-1, E, 1, EF, 1 );
         cpttrf(N, DF, EF, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = CLANHT( '1', N, D, E )

      // Compute the reciprocal of the condition number of A.

      cptcon(N, DF, EF, ANORM, RCOND, RWORK, INFO );

      // Compute the solution vectors X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      cpttrs('Lower', N, NRHS, DF, EF, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      cptrfs('Lower', N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of CPTSVX

      }
