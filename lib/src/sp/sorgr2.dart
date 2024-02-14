      void sorgr2(final int M, final int N, final int K, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, M, N;
      double               A( LDA, * ), TAU( * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, II, J, L;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( K < 0 || K > M ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SORGR2', -INFO );
         return;
      }

      // Quick return if possible

      if (M <= 0) return;

      if ( K < M ) {

         // Initialise rows 1:m-k to rows of the unit matrix

         for (J = 1; J <= N; J++) { // 20
            for (L = 1; L <= M - K; L++) { // 10
               A[L][J] = ZERO;
            } // 10
            if (J > N-M && J <= N-K) A( M-N+J, J ) = ONE;
         } // 20
      }

      for (I = 1; I <= K; I++) { // 40
         II = M - K + I;

         // Apply H(i) to A(1:m-k+i,1:n-k+i) from the right

         A[II][N-M+II] = ONE;
         slarf('Right', II-1, N-M+II, A( II, 1 ), LDA, TAU( I ), A, LDA, WORK );
         sscal(N-M+II-1, -TAU( I ), A( II, 1 ), LDA );
         A[II][N-M+II] = ONE - TAU( I );

         // Set A(m-k+i,n-k+i+1:n) to zero

         for (L = N - M + II + 1; L <= N; L++) { // 30
            A[II][L] = ZERO;
         } // 30
      } // 40
      }
