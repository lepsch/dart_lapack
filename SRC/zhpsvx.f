      SUBROUTINE ZHPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )
      COMPLEX*16         AFP( * ), AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOFACT;
      DOUBLE PRECISION   ANORM
*     ..
*     .. External Functions ..
      bool               LSAME;
      DOUBLE PRECISION   DLAMCH, ZLANHP
      EXTERNAL           LSAME, DLAMCH, ZLANHP
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZHPCON, ZHPRFS, ZHPTRF, ZHPTRS, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
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
         CALL XERBLA( 'ZHPSVX', -INFO )
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the factorization A = U*D*U**H or A = L*D*L**H.
*
         CALL ZCOPY( N*( N+1 ) / 2, AP, 1, AFP, 1 )
         CALL ZHPTRF( UPLO, N, AFP, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.GT.0 )THEN
            RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A.
*
      ANORM = ZLANHP( 'I', UPLO, N, AP, RWORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL ZHPCON( UPLO, N, AFP, IPIV, ANORM, RCOND, WORK, INFO )
*
*     Compute the solution vectors X.
*
      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZHPTRS( UPLO, N, NRHS, AFP, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
*     End of ZHPSVX
*
      END
