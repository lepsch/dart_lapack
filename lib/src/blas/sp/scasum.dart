      double scasum(final int N, final int CX, final int INCX,) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,N;
      Complex CX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double STEMP;
      int     I,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS,AIMAG,REAL
      // ..
      SCASUM = 0.0;
      STEMP = 0.0;
      if (N <= 0 || INCX <= 0) return;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            STEMP = STEMP + ABS(double(CX(I))) + ABS(AIMAG(CX(I)));
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
            STEMP = STEMP + ABS(double(CX(I))) + ABS(AIMAG(CX(I)));
         }
      }
      SCASUM = STEMP;
      }
