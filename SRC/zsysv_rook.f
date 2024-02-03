      SUBROUTINE ZSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSYTRF_ROOK, ZSYTRS_ROOK
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
         INFO = -8
      } else if ( LWORK.LT.1 .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            LWKOPT = 1
         } else {
            zsytrf_rook(UPLO, N, A, LDA, IPIV, WORK, -1, INFO );
            LWKOPT = INT( DBLE( WORK( 1 ) ) )
         }
         WORK( 1 ) = LWKOPT
      }

      if ( INFO.NE.0 ) {
         xerbla('ZSYSV_ROOK ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = U*D*U**T or A = L*D*L**T.

      zsytrf_rook(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         // Solve with TRS_ROOK ( Use Level 2 BLAS)

         zsytrs_rook(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO );

      }

      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZSYSV_ROOK

      }
