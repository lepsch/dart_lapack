      SUBROUTINE SSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS;
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, NOFACT;
      int                LWKMIN, LWKOPT, NB;
      REAL               ANORM
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANSY, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLACPY, SSYCON, SSYRFS, SSYTRF, SSYTRS, XERBLA
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
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 3*N )
      IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
*
      IF( INFO.EQ.0 ) THEN
         LWKOPT = LWKMIN
         IF( NOFACT ) THEN
            NB = ILAENV( 1, 'SSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = MAX( LWKOPT, N*NB )
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYSVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the factorization A = U*D*U**T or A = L*D*L**T.
*
         CALL SLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL SSYTRF( UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO )
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
      ANORM = SLANSY( 'I', UPLO, N, A, LDA, WORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL SSYCON( UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, IWORK, INFO )
*
*     Compute the solution vectors X.
*
      CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL SSYTRS( UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL SSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
      RETURN
*
*     End of SSYSVX
*
      END
