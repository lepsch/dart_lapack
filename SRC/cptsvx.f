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
      PARAMETER          ( ZERO = 0.0E+0 )
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
      IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPTSVX', -INFO )
         RETURN
      END IF

      IF( NOFACT ) THEN

         // Compute the L*D*L**H (or U**H*D*U) factorization of A.

         CALL SCOPY( N, D, 1, DF, 1 )
         IF( N.GT.1 ) CALL CCOPY( N-1, E, 1, EF, 1 )
         CALL CPTTRF( N, DF, EF, INFO )

         // Return if INFO is non-zero.

         IF( INFO.GT.0 )THEN
            RCOND = ZERO
            RETURN
         END IF
      END IF

      // Compute the norm of the matrix A.

      ANORM = CLANHT( '1', N, D, E )

      // Compute the reciprocal of the condition number of A.

      CALL CPTCON( N, DF, EF, ANORM, RCOND, RWORK, INFO )

      // Compute the solution vectors X.

      CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL CPTTRS( 'Lower', N, NRHS, DF, EF, X, LDX, INFO )

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      CALL CPTRFS( 'Lower', N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of CPTSVX

      }
