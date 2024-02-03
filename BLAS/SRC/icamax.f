      int     FUNCTION ICAMAX(N,CX,INCX);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL SMAX
      int     I,IX;
      // ..
      // .. External Functions ..
      REAL SCABS1
      // EXTERNAL SCABS1
      // ..
      ICAMAX = 0
      if (N.LT.1 || INCX.LE.0) RETURN;
      ICAMAX = 1
      if (N == 1) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = SCABS1(CX(1))
         for (I = 2; I <= N; I++) {
            if (SCABS1(CX(I)).GT.SMAX) {
               ICAMAX = I
               SMAX = SCABS1(CX(I))
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1
         SMAX = SCABS1(CX(1))
         IX = IX + INCX
         for (I = 2; I <= N; I++) {
            if (SCABS1(CX(IX)).GT.SMAX) {
               ICAMAX = I
               SMAX = SCABS1(CX(IX))
            }
            IX = IX + INCX
         }
      }
      RETURN

      // End of ICAMAX

      }
