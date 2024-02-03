      double           FUNCTION DSDOT(N,SX,INCX,SY,INCY);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*)
      // ..

*  Authors:
*  ========
*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)

*  =====================================================================

      // .. Local Scalars ..
      int     I,KX,KY,NS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      DSDOT = 0.0D0
      if (N.LE.0) RETURN;
      if (INCX == INCY .AND. INCX.GT.0) {

      // Code for equal, positive, non-unit increments.

         NS = N*INCX
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         }
      } else {

      // Code for unequal or nonpositive increments.

         KX = 1
         KY = 1
         if (INCX.LT.0) KX = 1 + (1-N)*INCX;
         if (INCY.LT.0) KY = 1 + (1-N)*INCY;
         for (I = 1; I <= N; I++) {
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         }
      }
      RETURN

      // End of DSDOT

      }
