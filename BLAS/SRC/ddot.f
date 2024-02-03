      double           FUNCTION DDOT(N,DX,INCX,DY,INCY);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      double           DX(*),DY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double           DTEMP;
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      DDOT = 0.0;
      DTEMP = 0.0;
      if (N <= 0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1


         // clean-up loop

         M = MOD(N,5);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               DTEMP = DTEMP + DX(I)*DY(I);
            }
            if (N < 5) {
               DDOT=DTEMP;
            return;
            }
         }
         MP1 = M + 1;
         DO I = MP1,N,5;
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4);
         }
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            DTEMP = DTEMP + DX(IX)*DY(IY);
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      DDOT = DTEMP;
      return;

      // End of DDOT

      }
