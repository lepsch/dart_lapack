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
      IF (INCX.EQ.1) THEN

         // code for increment equal to 1

         SMAX = SCABS1(CX(1))
         DO I = 2,N
            IF (SCABS1(CX(I)).GT.SMAX) THEN
               ICAMAX = I
               SMAX = SCABS1(CX(I))
            END IF
         END DO
      ELSE

         // code for increment not equal to 1

         IX = 1
         SMAX = SCABS1(CX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (SCABS1(CX(IX)).GT.SMAX) THEN
               ICAMAX = I
               SMAX = SCABS1(CX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN

      // End of ICAMAX

      }
