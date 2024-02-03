      int     FUNCTION IZMAX1( N, ZX, INCX );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         ZX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double             DMAX;
      int                I, IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      IZMAX1 = 0;
      if (N < 1 || INCX <= 0) RETURN;
      IZMAX1 = 1;
      if (N == 1) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         DMAX = ABS(ZX(1));
         for (I = 2; I <= N; I++) {
            if (ABS(ZX(I)) > DMAX) {
               IZMAX1 = I;
               DMAX = ABS(ZX(I));
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         DMAX = ABS(ZX(1));
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if (ABS(ZX(IX)) > DMAX) {
               IZMAX1 = I;
               DMAX = ABS(ZX(IX));
            }
            IX = IX + INCX;
         }
      }
      return;

      // End of IZMAX1

      }
