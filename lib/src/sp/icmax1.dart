// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      int icmax1(final int N, final int CX, final int INCX,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      Complex            CX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double               SMAX;
      int                I, IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      ICMAX1 = 0;
      if (N < 1 || INCX <= 0) return;
      ICMAX1 = 1;
      if (N == 1) return;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = (CX(1)).abs();
         for (I = 2; I <= N; I++) {
            if ((CX(I)).abs() > SMAX) {
               ICMAX1 = I;
               SMAX = (CX(I)).abs();
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         SMAX = (CX(1)).abs();
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if ((CX(IX)).abs() > SMAX) {
               ICMAX1 = I;
               SMAX = (CX(IX)).abs();
            }
            IX = IX + INCX;
         }
      }
      }
