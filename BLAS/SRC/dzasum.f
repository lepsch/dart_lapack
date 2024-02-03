      double dzasum(N,ZX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      Complex ZX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double           STEMP;
      int     I,NINCX;
      // ..
      // .. External Functions ..
      double           DCABS1;
      // EXTERNAL DCABS1
      // ..
      DZASUM = 0.0;
      STEMP = 0.0;
      if (N <= 0 || INCX <= 0) return;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            STEMP = STEMP + DCABS1(ZX(I));
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         DO I = 1,NINCX,INCX;
            STEMP = STEMP + DCABS1(ZX(I));
         }
      }
      DZASUM = STEMP;
      return;
      }
