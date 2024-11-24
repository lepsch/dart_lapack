// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      int ilaclr(final int M, final int N, final int A, final int LDA,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                M, N, LDA;
      Complex            A( LDA, * );
      // ..

      Complex          ZERO;
      const     ZERO = (0.0, 0.0) ;
      int     I, J;

      // Quick test for the common case where one corner is non-zero.
      if ( M == 0 ) {
         ILACLR = M;
      } else if ( A(M, 1) != ZERO || A(M, N) != ZERO ) {
         ILACLR = M;
      } else {
      // Scan up each column tracking the last zero row seen.
         ILACLR = 0;
         for (J = 1; J <= N; J++) {
            I=M;
            while ((A(max(I,1),J) == ZERO) && (I >= 1)) {
               I=I-1;
            }
            ILACLR = max( ILACLR, I );
         }
      }
      }
