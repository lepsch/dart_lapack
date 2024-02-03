      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      double           DX(*),DY(*);
      // ..

*  =====================================================================

      // .. Local Scalars ..
      double           DTEMP;
      int     I,IX,IY,M,MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      if (N.LE.0) RETURN;
      if (INCX.EQ.1 .AND. INCY.EQ.1) {

        // code for both increments equal to 1


        // clean-up loop

         M = MOD(N,3)
         if (M.NE.0) {
            for (I = 1; I <= M; I++) {
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            }
            if (N.LT.3) RETURN;
         }
         MP1 = M + 1
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         }
      } else {

        // code for unequal increments or equal increments not equal
          // to 1

         IX = 1
         IY = 1
         if (INCX.LT.0) IX = (-N+1)*INCX + 1;
         if (INCY.LT.0) IY = (-N+1)*INCY + 1;
         for (I = 1; I <= N; I++) {
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         }
      }
      RETURN

      // End of DSWAP

      }
