      void saxpy(N,SA,SX,INCX,SY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double SA;
      int     INCX,INCY,N;
      double SX(*),SY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      if (N <= 0) return;
      if (SA == 0.0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1


         // clean-up loop

         M = (N % 4);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               SY[I] = SY(I) + SA*SX(I);
            }
         }
         if (N < 4) return;
         MP1 = M + 1;
         for (I = MP1; I <= N; I += 4) {
            SY[I] = SY(I) + SA*SX(I);
            SY[I+1] = SY(I+1) + SA*SX(I+1);
            SY[I+2] = SY(I+2) + SA*SX(I+2);
            SY[I+3] = SY(I+3) + SA*SX(I+3);
         }
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
          SY[IY] = SY(IY) + SA*SX(IX);
          IX = IX + INCX;
          IY = IY + INCY;
         }
      }
      return;
      }
