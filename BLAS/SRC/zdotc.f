      COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)

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
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG
      // ..
      ZTEMP = (0.0d0,0.0d0)
      ZDOTC = (0.0d0,0.0d0)
      if (N.LE.0) RETURN;
      if (INCX == 1 .AND. INCY == 1) {

         // code for both increments equal to 1

         for (I = 1; I <= N; I++) {
            ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
         }
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1
         IY = 1
         if (INCX.LT.0) IX = (-N+1)*INCX + 1;
         if (INCY.LT.0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         }
      }
      ZDOTC = ZTEMP
      RETURN

      // End of ZDOTC

      }
