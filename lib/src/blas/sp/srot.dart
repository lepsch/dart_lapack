// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void srot(final int N, final int SX, final int INCX, final int SY, final int INCY, final int C, final int S,) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double C,S;
      int     INCX,INCY,N;
      double SX(*),SY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double STEMP;
      int     I,IX,IY;
      // ..
      if (N <= 0) return;
      if (INCX == 1 && INCY == 1) {

        // code for both increments equal to 1

         for (I = 1; I <= N; I++) {
            STEMP = C*SX(I) + S*SY(I);
            SY[I] = C*SY(I) - S*SX(I);
            SX[I] = STEMP;
         }
      } else {

        // code for unequal increments or equal increments not equal
        //   to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            STEMP = C*SX(IX) + S*SY(IY);
            SY[IY] = C*SY(IY) - S*SX(IX);
            SX[IX] = STEMP;
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      }
