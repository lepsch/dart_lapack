      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)

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
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      CTEMP = (0.0,0.0)
      CDOTC = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN

         // code for both increments equal to 1

         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(I))*CY(I)
         END DO
      ELSE

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CDOTC = CTEMP
      RETURN

      // End of CDOTC

      END
