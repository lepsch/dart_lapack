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
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ICAMAX = 1
      IF (N.EQ.1) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1

         SMAX = SCABS1(CX(1))
         DO I = 2,N
            if (SCABS1(CX(I)).GT.SMAX) {
               ICAMAX = I
               SMAX = SCABS1(CX(I))
            }
         END DO
      } else {

         // code for increment not equal to 1

         IX = 1
         SMAX = SCABS1(CX(1))
         IX = IX + INCX
         DO I = 2,N
            if (SCABS1(CX(IX)).GT.SMAX) {
               ICAMAX = I
               SMAX = SCABS1(CX(IX))
            }
            IX = IX + INCX
         END DO
      }
      RETURN

      // End of ICAMAX

      }
