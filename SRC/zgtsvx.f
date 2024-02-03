      SUBROUTINE ZGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             FACT, TRANS;
      int                INFO, LDB, LDX, N, NRHS
      double             RCOND;
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      double             BERR( * ), FERR( * ), RWORK( * );
      COMPLEX*16         B( LDB, * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOFACT, NOTRAN;
      String             NORM;
      double             ANORM;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGT;
      EXTERNAL           LSAME, DLAMCH, ZLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGTCON, ZGTRFS, ZGTTRF, ZGTTRS, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
         CALL XERBLA( 'ZGTSVX', -INFO )
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the LU factorization of A.
*
         CALL ZCOPY( N, D, 1, DF, 1 )
         IF( N.GT.1 ) THEN
            CALL ZCOPY( N-1, DL, 1, DLF, 1 )
            CALL ZCOPY( N-1, DU, 1, DUF, 1 )
         END IF
         CALL ZGTTRF( N, DLF, DF, DUF, DU2, IPIV, INFO )
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
      ANORM = ZLANGT( NORM, N, DL, D, DU )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL ZGTCON( NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, INFO )
*
*     Compute the solution vectors X.
*
      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZGTTRS( TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
*     End of ZGTSVX
*
      END
