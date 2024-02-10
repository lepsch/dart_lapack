      Complex cdotu(N,CX,INCX,CY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,INCY,N;
      Complex CX(*),CY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      Complex CTEMP;
      int     I,IX,IY;
      // ..
      CTEMP = (0.0,0.0);
      CDOTU = (0.0,0.0);
      if (N <= 0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1

         for (I = 1; I <= N; I++) {
            CTEMP = CTEMP + CX(I)*CY(I);
         }
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            CTEMP = CTEMP + CX(IX)*CY(IY);
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      CDOTU = CTEMP;
      }
