      REAL sdsdot(N,SB,SX,INCX,SY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL SB;
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*);
      // .. Local Scalars ..
      double           DSDOT;
      int     I,KX,KY,NS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      DSDOT = SB;
      if (N <= 0) {
         SDSDOT = REAL(DSDOT);
         return;
      }
      if (INCX == INCY && INCX > 0) {

      // Code for equal and positive increments.

         NS = N*INCX;
         DO I = 1,NS,INCX;
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I));
         }
      } else {

      // Code for unequal or nonpositive increments.

         KX = 1;
         KY = 1;
         if (INCX < 0) KX = 1 + (1-N)*INCX;
         if (INCY < 0) KY = 1 + (1-N)*INCY;
         for (I = 1; I <= N; I++) {
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY));
            KX = KX + INCX;
            KY = KY + INCY;
         }
      }
      SDSDOT = REAL(DSDOT);
      return;

      // End of SDSDOT

      }
