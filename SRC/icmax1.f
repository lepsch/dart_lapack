      int     FUNCTION ICMAX1( N, CX, INCX );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            CX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL               SMAX
      int                I, IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      ICMAX1 = 0
      if (N < 1 || INCX.LE.0) RETURN;
      ICMAX1 = 1
      if (N == 1) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = ABS(CX(1))
         for (I = 2; I <= N; I++) {
            if (ABS(CX(I)).GT.SMAX) {
               ICMAX1 = I
               SMAX = ABS(CX(I))
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1
         SMAX = ABS(CX(1))
         IX = IX + INCX
         for (I = 2; I <= N; I++) {
            if (ABS(CX(IX)).GT.SMAX) {
               ICMAX1 = I
               SMAX = ABS(CX(IX))
            }
            IX = IX + INCX
         }
      }
      RETURN

      // End of ICMAX1

      }
