      SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*),CY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      COMPLEX CTEMP
      int     I,IX,IY;
      // ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN

        // code for both increments equal to 1
         DO I = 1,N
            CTEMP = CX(I)
            CX(I) = CY(I)
            CY(I) = CTEMP
         END DO
      } else {

        // code for unequal increments or equal increments not equal
         t // o 1

         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CX(IX)
            CX(IX) = CY(IY)
            CY(IY) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN

      // End of CSWAP

      }
