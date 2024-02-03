      SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*),CY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      COMPLEX CTEMP
      int     I,IX,IY;
      // ..
      if (N.LE.0) RETURN;
      if (INCX == 1 && INCY == 1) {

        // code for both increments equal to 1
         for (I = 1; I <= N; I++) {
            CTEMP = CX(I)
            CX(I) = CY(I)
            CY(I) = CTEMP
         }
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1
         IY = 1
         if (INCX.LT.0) IX = (-N+1)*INCX + 1;
         if (INCY.LT.0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            CTEMP = CX(IX)
            CX(IX) = CY(IY)
            CY(IY) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         }
      }
      RETURN

      // End of CSWAP

      }
