      int     FUNCTION IZAMAX(N,ZX,INCX);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 ZX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      double           DMAX;
      int     I,IX;
      // ..
      // .. External Functions ..
      double           DCABS1;
      // EXTERNAL DCABS1
      // ..
      IZAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN

         // code for increment equal to 1

         DMAX = DCABS1(ZX(1))
         DO I = 2,N
            IF (DCABS1(ZX(I)).GT.DMAX) THEN
               IZAMAX = I
               DMAX = DCABS1(ZX(I))
            END IF
         END DO
      ELSE

         // code for increment not equal to 1

         IX = 1
         DMAX = DCABS1(ZX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DCABS1(ZX(IX)).GT.DMAX) THEN
               IZAMAX = I
               DMAX = DCABS1(ZX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN

      // End of IZAMAX

      END
