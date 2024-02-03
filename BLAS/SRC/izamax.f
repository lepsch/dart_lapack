      int     FUNCTION IZAMAX(N,ZX,INCX);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 ZX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double           DMAX;
      int     I,IX;
      // ..
      // .. External Functions ..
      double           DCABS1;
      // EXTERNAL DCABS1
      // ..
      IZAMAX = 0;
      if (N < 1 || INCX <= 0) return;
      IZAMAX = 1;
      if (N == 1) return;
      if (INCX == 1) {

         // code for increment equal to 1

         DMAX = DCABS1(ZX(1));
         for (I = 2; I <= N; I++) {
            if (DCABS1(ZX(I)) > DMAX) {
               IZAMAX = I;
               DMAX = DCABS1(ZX(I));
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         DMAX = DCABS1(ZX(1));
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if (DCABS1(ZX(IX)) > DMAX) {
               IZAMAX = I;
               DMAX = DCABS1(ZX(IX));
            }
            IX = IX + INCX;
         }
      }
      return;

      // End of IZAMAX

      }
