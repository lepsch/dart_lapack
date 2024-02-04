      int isamax(N,SX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double SX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double SMAX;
      int     I,IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      ISAMAX = 0;
      if (N < 1 || INCX <= 0) return;
      ISAMAX = 1;
      if (N == 1) return;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = (SX(1)).abs();
         for (I = 2; I <= N; I++) {
            if ((SX(I)).abs() > SMAX) {
               ISAMAX = I;
               SMAX = (SX(I)).abs();
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         SMAX = (SX(1)).abs();
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if ((SX(IX)).abs() > SMAX) {
               ISAMAX = I;
               SMAX = (SX(IX)).abs();
            }
            IX = IX + INCX;
         }
      }
      return;
      }
