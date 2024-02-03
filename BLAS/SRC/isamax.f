      int     FUNCTION ISAMAX(N,SX,INCX);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL SMAX
      int     I,IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      ISAMAX = 0
      if (N.LT.1 .OR. INCX.LE.0) RETURN;
      ISAMAX = 1
      if (N == 1) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = ABS(SX(1))
         for (I = 2; I <= N; I++) {
            if (ABS(SX(I)).GT.SMAX) {
               ISAMAX = I
               SMAX = ABS(SX(I))
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1
         SMAX = ABS(SX(1))
         IX = IX + INCX
         for (I = 2; I <= N; I++) {
            if (ABS(SX(IX)).GT.SMAX) {
               ISAMAX = I
               SMAX = ABS(SX(IX))
            }
            IX = IX + INCX
         }
      }
      RETURN

      // End of ISAMAX

      }
