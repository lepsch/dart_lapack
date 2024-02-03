      void daxpy(N,DA,DX,INCX,DY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           DA;
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      double           DX(*),DY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      if (N <= 0) return;
      if (DA == 0.0) return;
      if (INCX == 1 && INCY == 1) {

         // code for both increments equal to 1


         // clean-up loop

         M = (N % 4);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               DY(I) = DY(I) + DA*DX(I);
            }
         }
         if (N < 4) return;
         MP1 = M + 1;
         for (I = MP1; I <= N; I += 4) {
            DY(I) = DY(I) + DA*DX(I);
            DY(I+1) = DY(I+1) + DA*DX(I+1);
            DY(I+2) = DY(I+2) + DA*DX(I+2);
            DY(I+3) = DY(I+3) + DA*DX(I+3);
         }
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1;
         IY = 1;
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
          DY(IY) = DY(IY) + DA*DX(IX);
          IX = IX + INCX;
          IY = IY + INCY;
         }
      }
      return;
      }
