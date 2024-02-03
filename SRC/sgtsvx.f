      SUBROUTINE SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             FACT, TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               NOFACT, NOTRAN;
      String             NORM;
      REAL               ANORM
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGT
      EXTERNAL           LSAME, SLAMCH, SLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGTCON, SGTRFS, SGTTRF, SGTTRS, SLACPY, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGTSVX', -INFO )
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the LU factorization of A.
*
         CALL SCOPY( N, D, 1, DF, 1 )
         IF( N.GT.1 ) THEN
            CALL SCOPY( N-1, DL, 1, DLF, 1 )
            CALL SCOPY( N-1, DU, 1, DUF, 1 )
         END IF
         CALL SGTTRF( N, DLF, DF, DUF, DU2, IPIV, INFO )
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
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = SLANGT( NORM, N, DL, D, DU )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL SGTCON( NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
*
*     Compute the solution vectors X.
*
      CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL SGTTRS( TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL SGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
*     End of SGTSVX
*
      END
