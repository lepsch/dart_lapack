      void zptsv(N, NRHS, D, E, final Matrix<double> B, final int LDB, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDB, N, NRHS;
      double             D( * );
      Complex         B( LDB, * ), E( * );
      // ..

// =====================================================================

      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZPTTRF, ZPTTRS
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
         xerbla('ZPTSV ', -INFO );
         return;
      }

      // Compute the L*D*L**H (or U**H*D*U) factorization of A.

      zpttrf(N, D, E, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         zpttrs('Lower', N, NRHS, D, E, B, LDB, INFO );
      }
      }
