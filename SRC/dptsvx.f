      SUBROUTINE DPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             FACT;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
*     ..
*     .. Array Arguments ..
      double             B( LDB, * ), BERR( * ), D( * ), DF( * ), E( * ), EF( * ), FERR( * ), WORK( * ), X( LDX, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOFACT;
      double             ANORM;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANST;
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DPTCON, DPTRFS, DPTTRF, DPTTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
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
         CALL XERBLA( 'DPTSVX', -INFO )
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the L*D*L**T (or U**T*D*U) factorization of A.
*
         CALL DCOPY( N, D, 1, DF, 1 )
         IF( N.GT.1 ) CALL DCOPY( N-1, E, 1, EF, 1 )
         CALL DPTTRF( N, DF, EF, INFO )
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
      ANORM = DLANST( '1', N, D, E )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL DPTCON( N, DF, EF, ANORM, RCOND, WORK, INFO )
*
*     Compute the solution vectors X.
*
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DPTTRS( N, NRHS, DF, EF, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
*     End of DPTSVX
*
      END
