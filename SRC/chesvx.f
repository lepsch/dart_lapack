      SUBROUTINE CHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               BERR( * ), FERR( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOFACT;
      int                LWKMIN, LWKOPT, NB;
      REAL               ANORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANHE, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, CLANHE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHECON, CHERFS, CHETRF, CHETRS, CLACPY, XERBLA
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
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N )
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
            NB = ILAENV( 1, 'CHETRF', UPLO, N, -1, -1, -1 )
            LWKOPT = MAX( LWKOPT, N*NB )
         END IF
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHESVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
         // Compute the factorization A = U*D*U**H or A = L*D*L**H.
*
         CALL CLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL CHETRF( UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO )
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
      ANORM = CLANHE( 'I', UPLO, N, A, LDA, RWORK )
*
      // Compute the reciprocal of the condition number of A.
*
      CALL CHECON( UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, INFO )
*
      // Compute the solution vectors X.
*
      CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL CHETRS( UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.
*
      CALL CHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
      // Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      RETURN
*
      // End of CHESVX
*
      END
