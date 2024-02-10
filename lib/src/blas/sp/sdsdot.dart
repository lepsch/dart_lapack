      double sdsdot(N,SB,SX,INCX,SY, final int INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double SB;
      int     INCX,INCY,N;
      double SX(*),SY(*);
      // .. Local Scalars ..
      double           DSDOT;
      int     I,KX,KY,NS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      DSDOT = SB;
      if (N <= 0) {
         SDSDOT = double(DSDOT);
         return;
      }
      if (INCX == INCY && INCX > 0) {

      // Code for equal and positive increments.

         NS = N*INCX;
         for (I = 1; INCX < 0 ? I >= NS : I <= NS; I += INCX) {
            DSDOT = DSDOT + (SX(I)).toDouble()*(SY(I)).toDouble();
         }
      } else {

      // Code for unequal or nonpositive increments.

         KX = 1;
         KY = 1;
         if (INCX < 0) KX = 1 + (1-N)*INCX;
         if (INCY < 0) KY = 1 + (1-N)*INCY;
         for (I = 1; I <= N; I++) {
            DSDOT = DSDOT + (SX(KX)).toDouble()*(SY(KY)).toDouble();
            KX = KX + INCX;
            KY = KY + INCY;
         }
      }
      SDSDOT = double(DSDOT);
      }
