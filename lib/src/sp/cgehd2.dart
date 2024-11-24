// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cgehd2(final int N, final int ILO, final int IHI, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, INFO, LDA, N;
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      int                I;
      Complex            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN

      // Test the input parameters

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -2;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('CGEHD2', -INFO );
         return;
      }

      for (I = ILO; I <= IHI - 1; I++) { // 10

         // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

         ALPHA = A( I+1, I );
         clarfg(IHI-I, ALPHA, A( min( I+2, N ), I ), 1, TAU( I ) );
         A[I+1][I] = ONE;

         // Apply H(i) to A(1:ihi,i+1:ihi) from the right

         clarf('Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), A( 1, I+1 ), LDA, WORK );

         // Apply H(i)**H to A(i+1:ihi,i+1:n) from the left

         clarf('Left', IHI-I, N-I, A( I+1, I ), 1, CONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK );

         A[I+1][I] = ALPHA;
      } // 10

      }
