      SUBROUTINE CSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK;
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CSYTRF_RK, CSYTRS_3
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
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9;
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -11;
      }

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            LWKOPT = 1;
         } else {
            csytrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, -1, INFO );
            LWKOPT = INT( WORK( 1 ) );
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('CSYSV_RK ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Compute the factorization A = U*D*U**T or A = L*D*L**T.

      csytrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO );

      if ( INFO == 0 ) {

         // Solve the system A*X = B with BLAS3 solver, overwriting B with X.

         csytrs_3(UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO );

      }

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);

      return;

      // End of CSYSV_RK

      }
