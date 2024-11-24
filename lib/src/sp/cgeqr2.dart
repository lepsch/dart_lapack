// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cgeqr2(final int M, final int N, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      int                I, K;
      Complex            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CGEQR2', -INFO );
         return;
      }

      K = min( M, N );

      for (I = 1; I <= K; I++) { // 10

         // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

         clarfg(M-I+1, A( I, I ), A( min( I+1, M ), I ), 1, TAU( I ) );
         if ( I < N ) {

            // Apply H(i)**H to A(i:m,i+1:n) from the left

            ALPHA = A( I, I );
            A[I][I] = ONE;
            clarf('Left', M-I+1, N-I, A( I, I ), 1, CONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK );
            A[I][I] = ALPHA;
         }
      } // 10
      }
