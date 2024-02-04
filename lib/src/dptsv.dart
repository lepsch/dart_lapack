      void dptsv(N, NRHS, D, E, B, LDB, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), E( * );
      // ..

// =====================================================================

      // .. External Subroutines ..
      // EXTERNAL DPTTRF, DPTTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

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
         xerbla('DPTSV ', -INFO );
         return;
      }

      // Compute the L*D*L**T (or U**T*D*U) factorization of A.

      dpttrf(N, D, E, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         dpttrs(N, NRHS, D, E, B, LDB, INFO );
      }
      return;
      }