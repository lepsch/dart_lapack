      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL SA
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      IF (N.LE.0) RETURN
      IF (SA.EQ.0.0) RETURN
      if (INCX.EQ.1 .AND. INCY.EQ.1) {

         // code for both increments equal to 1


         // clean-up loop

         M = MOD(N,4)
         if (M.NE.0) {
            for (I = 1; I <= M; I++) {
               SY(I) = SY(I) + SA*SX(I)
            END DO
         }
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            SY(I) = SY(I) + SA*SX(I)
            SY(I+1) = SY(I+1) + SA*SX(I+1)
            SY(I+2) = SY(I+2) + SA*SX(I+2)
            SY(I+3) = SY(I+3) + SA*SX(I+3)
         END DO
      } else {

         // code for unequal increments or equal increments
           // not equal to 1

         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         for (I = 1; I <= N; I++) {
          SY(IY) = SY(IY) + SA*SX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      }
      RETURN

      // End of SAXPY

      }
