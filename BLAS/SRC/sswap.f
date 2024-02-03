      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*),SY(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL STEMP
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      if (N.LE.0) RETURN;
      if (INCX == 1 && INCY == 1) {

        // code for both increments equal to 1


        // clean-up loop

         M = MOD(N,3)
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               STEMP = SX(I)
               SX(I) = SY(I)
               SY(I) = STEMP
            }
            if (N < 3) RETURN;
         }
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
         }
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1
         IY = 1
         if (INCX < 0) IX = (-N+1)*INCX + 1;
         if (INCY < 0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            STEMP = SX(IX)
            SX(IX) = SY(IY)
            SY(IY) = STEMP
            IX = IX + INCX
            IY = IY + INCY
         }
      }
      RETURN

      // End of SSWAP

      }
