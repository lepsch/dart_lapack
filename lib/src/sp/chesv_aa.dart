      void chesv_aa(final int UPLO, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Array<int> IPIV, final Matrix<double> B, final int LDB, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      int                IPIV( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKMIN, LWKOPT, LWKOPT_HETRF, LWKOPT_HETRS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHETRF_AA, CHETRS_AA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      LWKMIN = max( 1, 2*N, 3*N-2 );
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -10;
      }

      if ( INFO == 0 ) {
         chetrf_aa(UPLO, N, A, LDA, IPIV, WORK, -1, INFO );
         LWKOPT_HETRF = INT( WORK( 1 ) );
         chetrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO );
         LWKOPT_HETRS = INT( WORK( 1 ) );
         LWKOPT = max( LWKMIN, LWKOPT_HETRF, LWKOPT_HETRS );
         WORK[1] = SROUNDUP_LWORK( LWKOPT );
      }

      if ( INFO != 0 ) {
         xerbla('CHESV_AA ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Compute the factorization A = U**H*T*U or A = L*T*L**H.

      chetrf_aa(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         chetrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO );

      }

      WORK[1] = SROUNDUP_LWORK( LWKOPT );

      }
