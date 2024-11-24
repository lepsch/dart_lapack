// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void stpqrt(final int M, final int N, final int L, final int NB, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Matrix<double> T_, final int LDT, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final T = T_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INFO, LDA, LDB, LDT, N, M, L, NB;
      double A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      int        I, IB, LB, MB, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPQRT2, STPRFB, XERBLA

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
         xerbla('STPQRT', -INFO );
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

         stpqrt2(MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB, T(1, I ), LDT, IINFO );

      // Update by applying H^H to B(:,I+IB:N) from the left

         if ( I+IB <= N ) {
            stprfb('L', 'T', 'F', 'C', MB, N-I-IB+1, IB, LB, B( 1, I ), LDB, T( 1, I ), LDT, A( I, I+IB ), LDA, B( 1, I+IB ), LDB, WORK, IB );
         }
      }
      }
