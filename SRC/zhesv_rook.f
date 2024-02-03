      SUBROUTINE ZHESV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

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
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHETRF_ROOK, ZHETRS_ROOK
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
         INFO = -8
      } else if ( LWORK.LT.1 .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            LWKOPT = 1
         } else {
            NB = ILAENV( 1, 'ZHETRF_ROOK', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
         }
         WORK( 1 ) = LWKOPT
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHESV_ROOK ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = U*D*U**H or A = L*D*L**H.

      zhetrf_rook(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         // Solve with TRS ( Use Level BLAS 2)

         zhetrs_rook(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO );

      }

      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZHESV_ROOK

      }
