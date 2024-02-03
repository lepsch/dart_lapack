      SUBROUTINE CHESV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKMIN, LWKOPT, LWKOPT_HETRF, LWKOPT_HETRS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHETRF_AA, CHETRS_AA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N, 3*N-2 )
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      if ( INFO.EQ.0 ) {
         CALL CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
         LWKOPT_HETRF = INT( WORK( 1 ) )
         CALL CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO )
         LWKOPT_HETRS = INT( WORK( 1 ) )
         LWKOPT = MAX( LWKMIN, LWKOPT_HETRF, LWKOPT_HETRS )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CHESV_AA ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = U**H*T*U or A = L*T*L**H.

      CALL CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         CALL CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

      }

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of CHESV_AA

      }
