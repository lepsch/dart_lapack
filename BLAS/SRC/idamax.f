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
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1

         DMAX = DABS(DX(1))
         DO I = 2,N
            if (DABS(DX(I)).GT.DMAX) {
               IDAMAX = I
               DMAX = DABS(DX(I))
            }
         END DO
      } else {

         // code for increment not equal to 1

         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            if (DABS(DX(IX)).GT.DMAX) {
               IDAMAX = I
               DMAX = DABS(DX(IX))
            }
            IX = IX + INCX
         END DO
      }
      RETURN

      // End of IDAMAX

      }
