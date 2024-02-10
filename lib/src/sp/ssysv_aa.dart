      void ssysv_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      int                IPIV( * );
      double               A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKMIN, LWKOPT, LWKOPT_SYTRF, LWKOPT_SYTRS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      double               SROUNDUP_LWORK;
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SSYTRS_AA, SSYTRF_AA
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
         ssytrf_aa(UPLO, N, A, LDA, IPIV, WORK, -1, INFO );
         LWKOPT_SYTRF = INT( WORK( 1 ) );
         ssytrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO );
         LWKOPT_SYTRS = INT( WORK( 1 ) );
         LWKOPT = max( LWKMIN, LWKOPT_SYTRF, LWKOPT_SYTRS );
         WORK[1] = SROUNDUP_LWORK( LWKOPT );
      }

      if ( INFO != 0 ) {
         xerbla('SSYSV_AA', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Compute the factorization A = U**T*T*U or A = L*T*L**T.

      ssytrf_aa(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         ssytrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO );

      }

      WORK[1] = SROUNDUP_LWORK( LWKOPT );

      }
