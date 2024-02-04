      double dsdot(N,SX,INCX,SY,INCY) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*);
      // ..

// Authors:
// ========
// Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
// Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)

// =====================================================================

      // .. Local Scalars ..
      int     I,KX,KY,NS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      DSDOT = 0.0;
      if (N <= 0) return;
      if (INCX == INCY && INCX > 0) {

      // Code for equal, positive, non-unit increments.

         NS = N*INCX;
         for (I = 1; INCX < 0 ? I >= NS : I <= NS; I += INCX) {
            DSDOT = DSDOT + (SX(I)).toDouble()*(SY(I)).toDouble();
         }
      } else {

      // Code for unequal or nonpositive increments.

         KX = 1;
         KY = 1;
         if (INCX < 0) KX = 1 + (1-N)*INCX;
         if (INCY < 0) KY = 1 + (1-N)*INCY;
         for (I = 1; I <= N; I++) {
            DSDOT = DSDOT + (SX(KX)).toDouble()*(SY(KY)).toDouble();
            KX = KX + INCX;
            KY = KY + INCY;
         }
      }
      return;
      }
