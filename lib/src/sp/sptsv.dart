      void sptsv(final int N, final int NRHS, final int D, final int E, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final B = B_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDB, N, NRHS;
      double               B( LDB, * ), D( * ), E( * );
      // ..

// =====================================================================

      // .. External Subroutines ..
      // EXTERNAL SPTTRF, SPTTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SPTSV ', -INFO );
         return;
      }

      // Compute the L*D*L**T (or U**T*D*U) factorization of A.

      spttrf(N, D, E, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         spttrs(N, NRHS, D, E, B, LDB, INFO );
      }
      }
