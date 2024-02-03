      void zsysv_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT, LWKOPT_SYTRF, LWKOPT_SYTRS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSYTRF_AA, ZSYTRS_AA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < max(2*N, 3*N-2) && !LQUERY ) {
         INFO = -10;
      }

      if ( INFO == 0 ) {
         zsytrf_aa(UPLO, N, A, LDA, IPIV, WORK, -1, INFO );
         LWKOPT_SYTRF = INT( WORK(1) );
         zsytrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO );
         LWKOPT_SYTRS = INT( WORK(1) );
         LWKOPT = max( LWKOPT_SYTRF, LWKOPT_SYTRS );
         WORK( 1 ) = LWKOPT;
      }

      if ( INFO != 0 ) {
         xerbla('ZSYSV_AA ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Compute the factorization A = U**T*T*U or A = L*T*L**T.

      zsytrf_aa(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         zsytrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO );

      }

      WORK( 1 ) = LWKOPT;

      return;

      // End of ZSYSV_AA

      }
