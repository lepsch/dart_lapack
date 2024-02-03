      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      COMPLEX*16 ZTEMP
      int     I,IX,IY;
      // ..
      IF (N.LE.0) RETURN
      if (INCX.EQ.1 .AND. INCY.EQ.1) {

        // code for both increments equal to 1
         for (I = 1; I <= N; I++) {
            ZTEMP = ZX(I)
            ZX(I) = ZY(I)
            ZY(I) = ZTEMP
         END DO
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         for (I = 1; I <= N; I++) {
            ZTEMP = ZX(IX)
            ZX(IX) = ZY(IY)
            ZY(IY) = ZTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      }
      RETURN

      // End of ZSWAP

      }
