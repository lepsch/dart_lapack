      void sscal(N,SA,SX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double SA;
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double SX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,M,MP1,NINCX;
      // ..
      // .. Parameters ..
      double ONE;
      const     ONE=1.0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      if (N <= 0 || INCX <= 0 || SA == ONE) return;
      if (INCX == 1) {

         // code for increment equal to 1


         // clean-up loop

         M = (N % 5);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               SX[I] = SA*SX(I);
            }
            if (N < 5) return;
         }
         MP1 = M + 1;
         for (I = MP1; I <= N; I += 5) {
            SX[I] = SA*SX(I);
            SX[I+1] = SA*SX(I+1);
            SX[I+2] = SA*SX(I+2);
            SX[I+3] = SA*SX(I+3);
            SX[I+4] = SA*SX(I+4);
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
            SX[I] = SA*SX(I);
         }
      }
      return;
      }
