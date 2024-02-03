      SUBROUTINE CSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          FACT, UPLO
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      REAL               BERR( * ), FERR( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, NOFACT
      int                LWKOPT, NB
      REAL               ANORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      REAL               CLANSY, SLAMCH, SROUNDUP_LWORK
      EXTERNAL           ILAENV, LSAME, CLANSY, SLAMCH, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY, CSYCON, CSYRFS, CSYTRF, CSYTRS, XERBLA
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
      LQUERY = ( LWORK.EQ.-1 )
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
      ELSE IF( LWORK.LT.MAX( 1, 2*N ) .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
*
      IF( INFO.EQ.0 ) THEN
         LWKOPT = MAX( 1, 2*N )
         IF( NOFACT ) THEN
            NB = ILAENV( 1, 'CSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = MAX( LWKOPT, N*NB )
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYSVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the factorization A = U*D*U**T or A = L*D*L**T.
*
         CALL CLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL CSYTRF( UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO )
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
      ANORM = CLANSY( 'I', UPLO, N, A, LDA, RWORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL CSYCON( UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, INFO )
*
*     Compute the solution vectors X.
*
      CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL CSYTRS( UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL CSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
      RETURN
*
*     End of CSYSVX
*
      END
