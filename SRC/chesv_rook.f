      SUBROUTINE CHESV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INFO, LDA, LDB, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      int                LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, ILAENV, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CHETRF_ROOK, CHETRS_ROOK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
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
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'CHETRF_ROOK', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHESV_ROOK ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Compute the factorization A = U*D*U**H or A = L*D*L**H.
*
      CALL CHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
*        Solve with TRS ( Use Level BLAS 2)
*
         CALL CHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
      RETURN
*
*     End of CHESV_ROOK
*
      END
