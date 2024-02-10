      void zgeqrt(M, N, NB, final Matrix<double> A, final int LDA, final Matrix<double> T, final int LDT, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INFO, LDA, LDT, M, N, NB;
      Complex A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      int        I, IB, IINFO, K;
      bool       USE_RECURSIVE_QR;
      const    USE_RECURSIVE_QR= true ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEQRT2, ZGEQRT3, ZLARFB, XERBLA

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NB < 1 || ( NB > min(M,N) && min(M,N) > 0 ) ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDT < NB ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('ZGEQRT', -INFO );
         return;
      }

      // Quick return if possible

      K = min( M, N );
      if (K == 0) return;

      // Blocked loop of length K

      for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
         IB = min( K-I+1, NB );

      // Compute the QR factorization of the current block A(I:M,I:I+IB-1)

         if ( USE_RECURSIVE_QR ) {
            zgeqrt3(M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO );
         } else {
            zgeqrt2(M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO );
         }
         if ( I+IB <= N ) {

      // Update by applying H**H to A(I:M,I+IB:N) from the left

            zlarfb('L', 'C', 'F', 'C', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, T( 1, I ), LDT, A( I, I+IB ), LDA, WORK , N-I-IB+1 );
         }
      }
      }
