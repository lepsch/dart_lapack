      REAL FUNCTION SCASUM(N,CX,INCX);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      REAL STEMP;
      int     I,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS,AIMAG,REAL
      // ..
      SCASUM = 0.0;
      STEMP = 0.0;
      if (N <= 0 || INCX <= 0) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)));
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         DO I = 1,NINCX,INCX;
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)));
         }
      }
      SCASUM = STEMP;
      return;

      // End of SCASUM

      }
