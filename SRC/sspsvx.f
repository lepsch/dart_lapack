      SUBROUTINE SSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               AFP( * ), AP( * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
      // ..
*
*  =====================================================================
*
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
      REAL               SLAMCH, SLANSP
      // EXTERNAL LSAME, SLAMCH, SLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLACPY, SSPCON, SSPRFS, SSPTRF, SSPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSPSVX', -INFO )
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
         // Compute the factorization A = U*D*U**T or A = L*D*L**T.
*
         CALL SCOPY( N*( N+1 ) / 2, AP, 1, AFP, 1 )
         CALL SSPTRF( UPLO, N, AFP, IPIV, INFO )
*
         // Return if INFO is non-zero.
*
         IF( INFO.GT.0 )THEN
            RCOND = ZERO
            RETURN
         END IF
      END IF
*
      // Compute the norm of the matrix A.
*
      ANORM = SLANSP( 'I', UPLO, N, AP, WORK )
*
      // Compute the reciprocal of the condition number of A.
*
      CALL SSPCON( UPLO, N, AFP, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
*
      // Compute the solution vectors X.
*
      CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL SSPTRS( UPLO, N, NRHS, AFP, IPIV, X, LDX, INFO )
*
      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.
*
      CALL SSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
      // Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
      // End of SSPSVX
*
      END
