      double           FUNCTION DZASUM(N,ZX,INCX);

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
      double           STEMP;
      int     I,NINCX;
      // ..
      // .. External Functions ..
      double           DCABS1;
      // EXTERNAL DCABS1
      // ..
      DZASUM = 0.0d0
      STEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN

         // code for increment equal to 1

         DO I = 1,N
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      ELSE

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      END IF
      DZASUM = STEMP
      RETURN

      // End of DZASUM

      END
