      SUBROUTINE ZHESV_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LDB, LTB, LWORK, INFO;
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), TB( * ), WORK( * )
*     ..
*
*  =====================================================================
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE  = ( 1.0D+0, 0.0D+0 ) )
*
*     .. Local Scalars ..
      bool               UPPER, TQUERY, WQUERY;
      int                LWKOPT, LWKMIN;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHETRF_AA_2STAGE, ZHETRS_AA_2STAGE
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      WQUERY = ( LWORK.EQ.-1 )
      TQUERY = ( LTB.EQ.-1 )
      LWKMIN = MAX( 1, N )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LTB.LT.MAX( 1, 4*N ) .AND. .NOT.TQUERY ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.WQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
         CALL ZHETRF_AA_2STAGE( UPLO, N, A, LDA, TB, -1, IPIV, IPIV2, WORK, -1, INFO )
         LWKOPT = MAX( LWKMIN, INT( WORK( 1 ) ) )
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHESV_AA_2STAGE', -INFO )
         RETURN
      ELSE IF( WQUERY .OR. TQUERY ) THEN
         RETURN
      END IF
*
*     Compute the factorization A = U**H*T*U or A = L*T*L**H.
*
      CALL ZHETRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL ZHETRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO )
*
      END IF
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZHESV_AA_2STAGE
*
      END
