      int     FUNCTION IDAMAX(N,DX,INCX);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double           DX(*);
      // ..

*  =====================================================================

      // .. Local Scalars ..
      double           DMAX;
      int     I,IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS
      // ..
      IDAMAX = 0
      if (N.LT.1 || INCX.LE.0) RETURN;
      IDAMAX = 1
      if (N == 1) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         DMAX = DABS(DX(1))
         for (I = 2; I <= N; I++) {
            if (DABS(DX(I)).GT.DMAX) {
               IDAMAX = I
               DMAX = DABS(DX(I))
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         for (I = 2; I <= N; I++) {
            if (DABS(DX(IX)).GT.DMAX) {
               IDAMAX = I
               DMAX = DABS(DX(IX))
            }
            IX = IX + INCX
         }
      }
      RETURN

      // End of IDAMAX

      }
