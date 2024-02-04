      void stplqt(M, N, L, MB, A, LDA, B, LDB, T, LDT, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDB, LDT, N, M, L, MB;
      // ..
      // .. Array Arguments ..
      double A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, LB, NB, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPLQT2, STPRFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( L < 0 || (L > min(M,N) && min(M,N) >= 0)) {
         INFO = -3;
      } else if ( MB < 1 || (MB > M && M > 0)) {
         INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -6;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -8;
      } else if ( LDT < MB ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('STPLQT', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      for (I = 1; MB < 0 ? I >= M : I <= M; I += MB) {

      // Compute the QR factorization of the current block

         IB = min( M-I+1, MB );
         NB = min( N-L+I+IB-1, N );
         if ( I >= L ) {
            LB = 0;
         } else {
            LB = NB-N+L-I+1;
         }

         stplqt2(IB, NB, LB, A(I,I), LDA, B( I, 1 ), LDB, T(1, I ), LDT, IINFO );

      // Update by applying H**T to B(I+IB:M,:) from the right

         if ( I+IB <= M ) {
            stprfb('R', 'N', 'F', 'R', M-I-IB+1, NB, IB, LB, B( I, 1 ), LDB, T( 1, I ), LDT, A( I+IB, I ), LDA, B( I+IB, 1 ), LDB, WORK, M-I-IB+1);
         }
      }
      return;
      }