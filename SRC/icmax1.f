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
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ICMAX1 = 1
      IF (N.EQ.1) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1

         SMAX = ABS(CX(1))
         DO I = 2,N
            if (ABS(CX(I)).GT.SMAX) {
               ICMAX1 = I
               SMAX = ABS(CX(I))
            }
         END DO
      } else {

         // code for increment not equal to 1

         IX = 1
         SMAX = ABS(CX(1))
         IX = IX + INCX
         DO I = 2,N
            if (ABS(CX(IX)).GT.SMAX) {
               ICMAX1 = I
               SMAX = ABS(CX(IX))
            }
            IX = IX + INCX
         END DO
      }
      RETURN

      // End of ICMAX1

      }
