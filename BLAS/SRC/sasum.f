      REAL sasum(N,SX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      REAL STEMP;
      int     I,M,MP1,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS,MOD
      // ..
      SASUM = 0.0;
      STEMP = 0.0;
      if (N <= 0 || INCX <= 0) return;
      if (INCX == 1) {
         // code for increment equal to 1


         // clean-up loop

         M = (N % 6);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               STEMP = STEMP + (SX(I)).abs();
            }
            if (N < 6) {
               SASUM = STEMP;
               return;
            }
         }
         MP1 = M + 1;
         for (I = MP1; I <= N; I += 6) { //
            STEMP = STEMP + (SX(I)).abs() + (SX(I+1)).abs() + (SX(I+2)).abs() + (SX(I+3)).abs() + (SX(I+4)).abs() + (SX(I+5)).abs();
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) { //
            STEMP = STEMP + (SX(I)).abs();
         }
      }
      SASUM = STEMP;
      return;
      }
