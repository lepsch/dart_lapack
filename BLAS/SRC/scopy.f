      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
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
*  =====================================================================
*
*     .. Local Scalars ..
      int     I,IX,IY,M,MP1;
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               SY(I) = SX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,7
            SY(I) = SX(I)
            SY(I+1) = SX(I+1)
            SY(I+2) = SX(I+2)
            SY(I+3) = SX(I+3)
            SY(I+4) = SX(I+4)
            SY(I+5) = SX(I+5)
            SY(I+6) = SX(I+6)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            SY(IY) = SX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
*
*     End of SCOPY
*
      END
