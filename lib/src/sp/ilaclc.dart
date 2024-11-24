// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      int ilaclc(final int M, final int N, final int A, final int LDA,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                M, N, LDA;
      Complex            A( LDA, * );
      // ..

      Complex          ZERO;
      const     ZERO = (0.0, 0.0) ;
      int     I;

      // Quick test for the common case where one corner is non-zero.
      if ( N == 0 ) {
         ILACLC = N;
      } else if ( A(1, N) != ZERO || A(M, N) != ZERO ) {
         ILACLC = N;
      } else {
      // Now scan each column from the end, returning with the first non-zero.
         for (ILACLC = N; ILACLC >= 1; ILACLC--) {
            for (I = 1; I <= M; I++) {
               if( A(I, ILACLC) != ZERO ) return;
            }
         }
      }
      }
