      int     FUNCTION IZMAX1( N, ZX, INCX );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         ZX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      double             DMAX;
      int                I, IX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      IZMAX1 = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZMAX1 = 1
      IF (N.EQ.1) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1

         DMAX = ABS(ZX(1))
         DO I = 2,N
            if (ABS(ZX(I)).GT.DMAX) {
               IZMAX1 = I
               DMAX = ABS(ZX(I))
            }
         END DO
      } else {

         // code for increment not equal to 1

         IX = 1
         DMAX = ABS(ZX(1))
         IX = IX + INCX
         DO I = 2,N
            if (ABS(ZX(IX)).GT.DMAX) {
               IZMAX1 = I
               DMAX = ABS(ZX(IX))
            }
            IX = IX + INCX
         END DO
      }
      RETURN

      // End of IZMAX1

      }
