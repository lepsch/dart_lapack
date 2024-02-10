      void ctpqrt(M, N, L, NB, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> T, final int LDT, WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INFO, LDA, LDB, LDT, N, M, L, NB;
      Complex A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      int        I, IB, LB, MB, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTPQRT2, CTPRFB, XERBLA

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( L < 0 || (L > min(M,N) && min(M,N) >= 0)) {
         INFO = -3;
      } else if ( NB < 1 || (NB > N && N > 0)) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -8;
      } else if ( LDT < NB ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CTPQRT', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {

      // Compute the QR factorization of the current block

         IB = min( N-I+1, NB );
         MB = min( M-L+I+IB-1, M );
         if ( I >= L ) {
            LB = 0;
         } else {
            LB = MB-M+L-I+1;
         }

         ctpqrt2(MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB, T(1, I ), LDT, IINFO );

      // Update by applying H**H to B(:,I+IB:N) from the left

         if ( I+IB <= N ) {
            ctprfb('L', 'C', 'F', 'C', MB, N-I-IB+1, IB, LB, B( 1, I ), LDB, T( 1, I ), LDT, A( I, I+IB ), LDA, B( 1, I+IB ), LDB, WORK, IB );
         }
      }
      }
