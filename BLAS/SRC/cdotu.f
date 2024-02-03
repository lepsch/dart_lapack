      COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*),CY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      COMPLEX CTEMP;
      int     I,IX,IY;
      // ..
      CTEMP = (0.0,0.0);
      CDOTU = (0.0,0.0);
      if (N <= 0) RETURN;
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
      return;

      // End of CDOTU

      }
