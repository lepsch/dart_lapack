      SUBROUTINE CHESV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), E( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHETRF_RK, CHETRS_3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LWORK.LT.1 .AND. .NOT.LQUERY ) {
         INFO = -11
      }

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            LWKOPT = 1
         } else {
            chetrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, -1, INFO );
            LWKOPT = INT( WORK( 1 ) )
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         xerbla('CHESV_RK ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = U*D*U**T or A = L*D*L**T.

      chetrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO );

      if ( INFO == 0 ) {

         // Solve the system A*X = B with BLAS3 solver, overwriting B with X.

         chetrs_3(UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO );

      }

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

      RETURN

      // End of CHESV_RK

      }
