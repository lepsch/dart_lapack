      REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL SB
      int     INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     .. Local Scalars ..
      DOUBLE PRECISION DSDOT
      int     I,KX,KY,NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DBLE
*     ..
      DSDOT = SB
      IF (N.LE.0) THEN
         SDSDOT = REAL(DSDOT)
         RETURN
      END IF
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
*
*     Code for equal and positive increments.
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
      SDSDOT = REAL(DSDOT)
      RETURN
*
*     End of SDSDOT
*
      END
