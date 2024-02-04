      void srotm(N,SX,INCX,SY,INCY,SPARAM) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      double SPARAM(5),SX(*),SY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double SFLAG,SH11,SH12,SH21,SH22,TWO,W,Z,ZERO;
      int     I,KX,KY,NSTEPS;
      // ..
      // .. Data statements ..
      final (ZERO,TWO) = (0.0,2.0);
      // ..

      SFLAG = SPARAM(1);
      if (N <= 0 || (SFLAG+TWO == ZERO)) return;
      if (INCX == INCY && INCX > 0) {

         NSTEPS = N*INCX;
         if (SFLAG < ZERO) {
            SH11 = SPARAM(2);
            SH12 = SPARAM(4);
            SH21 = SPARAM(3);
            SH22 = SPARAM(5);
            for (I = 1; INCX < 0 ? I >= NSTEPS : I <= NSTEPS; I += INCX) {
               W = SX(I);
               Z = SY(I);
               SX[I] = W*SH11 + Z*SH12;
               SY[I] = W*SH21 + Z*SH22;
            }
         } else if (SFLAG == ZERO) {
            SH12 = SPARAM(4);
            SH21 = SPARAM(3);
            for (I = 1; INCX < 0 ? I >= NSTEPS : I <= NSTEPS; I += INCX) {
               W = SX(I);
               Z = SY(I);
               SX[I] = W + Z*SH12;
               SY[I] = W*SH21 + Z;
            }
         } else {
            SH11 = SPARAM(2);
            SH22 = SPARAM(5);
            for (I = 1; INCX < 0 ? I >= NSTEPS : I <= NSTEPS; I += INCX) {
               W = SX(I);
               Z = SY(I);
               SX[I] = W*SH11 + Z;
               SY[I] = -W + SH22*Z;
            }
         }
      } else {
         KX = 1;
         KY = 1;
         if (INCX < 0) KX = 1 + (1-N)*INCX;
         if (INCY < 0) KY = 1 + (1-N)*INCY;

         if (SFLAG < ZERO) {
            SH11 = SPARAM(2);
            SH12 = SPARAM(4);
            SH21 = SPARAM(3);
            SH22 = SPARAM(5);
            for (I = 1; I <= N; I++) {
               W = SX(KX);
               Z = SY(KY);
               SX[KX] = W*SH11 + Z*SH12;
               SY[KY] = W*SH21 + Z*SH22;
               KX = KX + INCX;
               KY = KY + INCY;
            }
         } else if (SFLAG == ZERO) {
            SH12 = SPARAM(4);
            SH21 = SPARAM(3);
            for (I = 1; I <= N; I++) {
               W = SX(KX);
               Z = SY(KY);
               SX[KX] = W + Z*SH12;
               SY[KY] = W*SH21 + Z;
               KX = KX + INCX;
               KY = KY + INCY;
            }
         } else {
             SH11 = SPARAM(2);
             SH22 = SPARAM(5);
             for (I = 1; I <= N; I++) {
                W = SX(KX);
                Z = SY(KY);
                SX[KX] = W*SH11 + Z;
                SY[KY] = -W + SH22*Z;
                KX = KX + INCX;
                KY = KY + INCY;
            }
         }
      }
      return;
      }
