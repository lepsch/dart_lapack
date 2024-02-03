      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*)
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      REAL STEMP
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
        // code for both increments equal to 1
*
*
        // clean-up loop
*
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               STEMP = SX(I)
               SX(I) = SY(I)
               SY(I) = STEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            STEMP = SX(I)
            SX(I) = SY(I)
            SY(I) = STEMP
            STEMP = SX(I+1)
            SX(I+1) = SY(I+1)
            SY(I+1) = STEMP
            STEMP = SX(I+2)
            SX(I+2) = SY(I+2)
            SY(I+2) = STEMP
         END DO
      ELSE
*
        // code for unequal increments or equal increments not equal
         t // o 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            STEMP = SX(IX)
            SX(IX) = SY(IY)
            SY(IY) = STEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
*
      // End of SSWAP
*
      END
