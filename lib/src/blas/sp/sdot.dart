// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      double sdot(final int N, final int SX, final int INCX, final int SY, final int INCY,) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,INCY,N;
      double SX(*),SY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double STEMP;
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      STEMP = 0.0;
      SDOT = 0.0;
      if (N <= 0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1


         // clean-up loop

         M = (N % 5);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               STEMP = STEMP + SX(I)*SY(I);
            }
            if (N < 5) {
               SDOT=STEMP;
            return;
            }
         }
         MP1 = M + 1;
         for (I = MP1; I <= N; I += 5) {
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4);
         }
      } else {

         // code for unequal increments or equal increments
         //   not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            STEMP = STEMP + SX(IX)*SY(IY);
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      SDOT = STEMP;
      }
