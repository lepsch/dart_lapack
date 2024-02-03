      SUBROUTINE CHESV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
*     ..
*     .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               LQUERY;
      int                LWKMIN, LWKOPT, LWKOPT_HETRF, LWKOPT_HETRS;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, ILAENV, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CHETRF_AA, CHETRS_AA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N, 3*N-2 )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         CALL CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
         LWKOPT_HETRF = INT( WORK( 1 ) )
         CALL CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO )
         LWKOPT_HETRS = INT( WORK( 1 ) )
         LWKOPT = MAX( LWKMIN, LWKOPT_HETRF, LWKOPT_HETRS )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHESV_AA ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Compute the factorization A = U**H*T*U or A = L*T*L**H.
*
      CALL CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
*
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      RETURN
*
*     End of CHESV_AA
*
      END
