      SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SPARAM(5),SX(*),SY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL SFLAG,SH11,SH12,SH21,SH22,TWO,W,Z,ZERO
      int     I,KX,KY,NSTEPS;
      // ..
      // .. Data statements ..
      DATA ZERO,TWO/0.E0,2.E0/
      // ..

      SFLAG = SPARAM(1)
      IF (N.LE.0 .OR. (SFLAG+TWO.EQ.ZERO)) RETURN
      IF (INCX.EQ.INCY.AND.INCX.GT.0) THEN

         NSTEPS = N*INCX
         IF (SFLAG.LT.ZERO) THEN
            SH11 = SPARAM(2)
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            SH22 = SPARAM(5)
            DO I = 1,NSTEPS,INCX
               W = SX(I)
               Z = SY(I)
               SX(I) = W*SH11 + Z*SH12
               SY(I) = W*SH21 + Z*SH22
            END DO
         ELSE IF (SFLAG.EQ.ZERO) THEN
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            DO I = 1,NSTEPS,INCX
               W = SX(I)
               Z = SY(I)
               SX(I) = W + Z*SH12
               SY(I) = W*SH21 + Z
            END DO
         ELSE
            SH11 = SPARAM(2)
            SH22 = SPARAM(5)
            DO I = 1,NSTEPS,INCX
               W = SX(I)
               Z = SY(I)
               SX(I) = W*SH11 + Z
               SY(I) = -W + SH22*Z
            END DO
         END IF
      ELSE
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY

         IF (SFLAG.LT.ZERO) THEN
            SH11 = SPARAM(2)
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            SH22 = SPARAM(5)
            DO I = 1,N
               W = SX(KX)
               Z = SY(KY)
               SX(KX) = W*SH11 + Z*SH12
               SY(KY) = W*SH21 + Z*SH22
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE IF (SFLAG.EQ.ZERO) THEN
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            DO I = 1,N
               W = SX(KX)
               Z = SY(KY)
               SX(KX) = W + Z*SH12
               SY(KY) = W*SH21 + Z
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE
             SH11 = SPARAM(2)
             SH22 = SPARAM(5)
             DO I = 1,N
                W = SX(KX)
                Z = SY(KY)
                SX(KX) = W*SH11 + Z
                SY(KY) = -W + SH22*Z
                KX = KX + INCX
                KY = KY + INCY
            END DO
         END IF
      END IF
      RETURN

      // End of SROTM

      END
