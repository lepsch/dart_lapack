      void chpsv(UPLO, N, NRHS, AP, final Array<int> IPIV, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      int                IPIV( * );
      Complex            AP( * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPTRF, CHPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CHPSV ', -INFO );
         return;
      }

      // Compute the factorization A = U*D*U**H or A = L*D*L**H.

      chptrf(UPLO, N, AP, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         chptrs(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO );

      }
      }
