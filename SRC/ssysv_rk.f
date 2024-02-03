      SUBROUTINE SSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), B( LDB, * ), E( * ), WORK( * )
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
      // EXTERNAL XERBLA, SSYTRF_RK, SSYTRS_3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
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

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            LWKOPT = 1
         } else {
            ssytrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, -1, INFO );
            LWKOPT = INT( WORK( 1 ) )
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         xerbla('SSYSV_RK ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = P*U*D*(U**T)*(P**T) or
      // A = P*U*D*(U**T)*(P**T).

      ssytrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO );

      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B with BLAS3 solver, overwriting B with X.

         ssytrs_3(UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO );

      }

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

      RETURN

      // End of SSYSV_RK

      }
