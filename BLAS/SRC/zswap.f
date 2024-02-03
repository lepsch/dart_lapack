      void zswap(N,ZX,INCX,ZY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      Complex ZX(*),ZY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      Complex ZTEMP;
      int     I,IX,IY;
      // ..
      if (N <= 0) return;
      if (INCX == 1 && INCY == 1) {

        // code for both increments equal to 1
         for (I = 1; I <= N; I++) {
            ZTEMP = ZX(I);
            ZX(I) = ZY(I);
            ZY(I) = ZTEMP;
         }
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            ZTEMP = ZX(IX);
            ZX(IX) = ZY(IY);
            ZY(IY) = ZTEMP;
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      return;

      // End of ZSWAP

      }
