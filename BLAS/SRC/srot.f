      SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL C,S;
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*);
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL STEMP;
      int     I,IX,IY;
      // ..
      if (N <= 0) RETURN;
      if (INCX == 1 && INCY == 1) {

        // code for both increments equal to 1

         for (I = 1; I <= N; I++) {
            STEMP = C*SX(I) + S*SY(I);
            SY(I) = C*SY(I) - S*SX(I);
            SX(I) = STEMP;
         }
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            STEMP = C*SX(IX) + S*SY(IY);
            SY(IY) = C*SY(IY) - S*SX(IX);
            SX(IX) = STEMP;
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      return;

      // End of SROT

      }
