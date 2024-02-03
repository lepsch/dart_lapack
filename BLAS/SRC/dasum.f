      double           FUNCTION DASUM(N,DX,INCX);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double           DX(*);
      // ..

*  =====================================================================

      // .. Local Scalars ..
      double           DTEMP;
      int     I,M,MP1,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS,MOD
      // ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      if (N <= 0 || INCX <= 0) RETURN;
      if (INCX == 1) {
         // code for increment equal to 1


         // clean-up loop

         M = MOD(N,6)
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               DTEMP = DTEMP + DABS(DX(I))
            }
            if (N < 6) {
               DASUM = DTEMP
               RETURN
            }
         }
         MP1 = M + 1
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         }
      }
      DASUM = DTEMP
      RETURN

      // End of DASUM

      }
