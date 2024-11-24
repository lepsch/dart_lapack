// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void csscal(final int N, final int SA, final int CX, final int INCX,) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double SA;
      int     INCX,N;
      Complex CX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,NINCX;
      // ..
      // .. Parameters ..
      double ONE;
      const     ONE=1.0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG,CMPLX,REAL
      // ..
      if (N <= 0 || INCX <= 0 || SA == ONE) return;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            CX[I] = CMPLX(SA*double(CX(I)),SA*AIMAG(CX(I)));
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
            CX[I] = CMPLX(SA*double(CX(I)),SA*AIMAG(CX(I)));
         }
      }
      }
