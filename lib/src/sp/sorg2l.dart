// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sorg2l(final int M, final int N, final int K, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> _WORK_, final Box<int> INFO,) {
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
      } else if ( N < 0 || N > M ) {
         INFO = -2;
      } else if ( K < 0 || K > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SORG2L', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) return;

      // Initialise columns 1:n-k to columns of the unit matrix

      for (J = 1; J <= N - K; J++) { // 20
         for (L = 1; L <= M; L++) { // 10
            A[L][J] = ZERO;
         } // 10
         A[M-N+J][J] = ONE;
      } // 20

      for (I = 1; I <= K; I++) { // 40
         II = N - K + I;

         // Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

         A[M-N+II][II] = ONE;
         slarf('Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, LDA, WORK );
         sscal(M-N+II-1, -TAU( I ), A( 1, II ), 1 );
         A[M-N+II][II] = ONE - TAU( I );

         // Set A(m-k+i+1:m,n-k+i) to zero

         for (L = M - N + II + 1; L <= M; L++) { // 30
            A[L][II] = ZERO;
         } // 30
      } // 40
      }
