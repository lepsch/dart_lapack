      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           C,S;
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      double           DX(*),DY(*);
      // ..

*  =====================================================================

      // .. Local Scalars ..
      double           DTEMP;
      int     I,IX,IY;
      // ..
      if (N <= 0) RETURN;
      if (INCX == 1 && INCY == 1) {

        // code for both increments equal to 1

         for (I = 1; I <= N; I++) {
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         }
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1
         IY = 1
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         }
      }
      RETURN

      // End of DROT

      }
