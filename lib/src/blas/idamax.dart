      int idamax(N,DX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double           DX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double           DMAX;
      int     I,IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS
      // ..
      IDAMAX = 0;
      if (N < 1 || INCX <= 0) return;
      IDAMAX = 1;
      if (N == 1) return;
      if (INCX == 1) {

         // code for increment equal to 1

         DMAX = (DX(1)).abs();
         for (I = 2; I <= N; I++) {
            if ((DX(I)).abs() > DMAX) {
               IDAMAX = I;
               DMAX = (DX(I)).abs();
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         DMAX = (DX(1)).abs();
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if ((DX(IX)).abs() > DMAX) {
               IDAMAX = I;
               DMAX = (DX(IX)).abs();
            }
            IX = IX + INCX;
         }
      }
      return;
      }
