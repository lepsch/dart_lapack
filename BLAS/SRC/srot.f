      SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL C,S
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL STEMP
      int     I,IX,IY;
      // ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN

        // code for both increments equal to 1

         DO I = 1,N
            STEMP = C*SX(I) + S*SY(I)
            SY(I) = C*SY(I) - S*SX(I)
            SX(I) = STEMP
         END DO
      ELSE

        // code for unequal increments or equal increments not equal
         t // o 1

         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            STEMP = C*SX(IX) + S*SY(IY)
            SY(IY) = C*SY(IY) - S*SX(IX)
            SX(IX) = STEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN

      // End of SROT

      }
