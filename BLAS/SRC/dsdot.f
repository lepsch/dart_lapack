      double           FUNCTION DSDOT(N,SX,INCX,SY,INCY);
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int     INCX,INCY,N;
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  Authors:
*  ========
*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
*
*  =====================================================================
*
*     .. Local Scalars ..
      int     I,KX,KY,NS;
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE
*     ..
      DSDOT = 0.0D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
*
*     Code for equal, positive, non-unit increments.
*
         NS = N*INCX
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         END DO
      ELSE
*
*     Code for unequal or nonpositive increments.
*
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
         DO I = 1,N
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         END DO
      END IF
      RETURN
*
*     End of DSDOT
*
      END
