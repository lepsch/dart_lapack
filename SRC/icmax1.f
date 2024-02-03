      int     FUNCTION ICMAX1( N, CX, INCX );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            CX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      REAL               SMAX;
      int                I, IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      ICMAX1 = 0;
      if (N < 1 || INCX <= 0) return;
      ICMAX1 = 1;
      if (N == 1) return;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = ABS(CX(1));
         for (I = 2; I <= N; I++) {
            if (ABS(CX(I)) > SMAX) {
               ICMAX1 = I;
               SMAX = ABS(CX(I));
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         SMAX = ABS(CX(1));
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if (ABS(CX(IX)) > SMAX) {
               ICMAX1 = I;
               SMAX = ABS(CX(IX));
            }
            IX = IX + INCX;
         }
      }
      return;

      // End of ICMAX1

      }
