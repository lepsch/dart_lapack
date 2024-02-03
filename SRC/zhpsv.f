      SUBROUTINE ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AP( * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHPTRF, ZHPTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('ZHPSV ', -INFO );
         RETURN
      }

      // Compute the factorization A = U*D*U**H or A = L*D*L**H.

      zhptrf(UPLO, N, AP, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         zhptrs(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO );

      }
      RETURN

      // End of ZHPSV

      }
