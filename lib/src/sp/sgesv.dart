      void sgesv(final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();
  final B = B_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, N, NRHS;
      int                IPIV( * );
      double               A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Subroutines ..
      // EXTERNAL SGETRF, SGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('SGESV ', -INFO );
         return;
      }

      // Compute the LU factorization of A.

      sgetrf(N, N, A, LDA, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         sgetrs('No transpose', N, NRHS, A, LDA, IPIV, B, LDB, INFO );
      }
      }
