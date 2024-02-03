      REAL sdot(N,SX,INCX,SY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      REAL STEMP;
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      STEMP = 0.0;
      SDOT = 0.0;
      if (N <= 0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1


         // clean-up loop

         M = MOD(N,5);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               STEMP = STEMP + SX(I)*SY(I);
            }
            if (N < 5) {
               SDOT=STEMP;
            return;
            }
         }
         MP1 = M + 1;
         DO I = MP1,N,5;
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4);
         }
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            STEMP = STEMP + SX(IX)*SY(IY);
            IX = IX + INCX;
            IY = IY + INCY;
         }
      }
      SDOT = STEMP;
      return;
      }
