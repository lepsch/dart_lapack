// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void caxpy(final int N, final int CA, final int CX, final int INCX, final int CY, final int INCY,) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      Complex CA;
      int     INCX,INCY,N;
      Complex CX(*),CY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,IX,IY;
      // ..
      // .. External Functions ..
      //- REAL SCABS1;
      // EXTERNAL SCABS1
      // ..
      if (N <= 0) return;
      if (SCABS1(CA) == 0.0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1

         for (I = 1; I <= N; I++) {
            CY[I] = CY(I) + CA*CX(I);
         }
      } else {

         // code for unequal increments or equal increments
         //   not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            CY[IY] = CY(IY) + CA*CX(IX);
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }

      }
