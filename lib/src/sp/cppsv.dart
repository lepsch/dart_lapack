      void cppsv(UPLO, N, NRHS, AP, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      Complex            AP( * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CPPTRF, CPPTRS, XERBLA
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
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CPPSV ', -INFO );
         return;
      }

      // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

      cpptrf(UPLO, N, AP, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         cpptrs(UPLO, N, NRHS, AP, B, LDB, INFO );

      }
      }
