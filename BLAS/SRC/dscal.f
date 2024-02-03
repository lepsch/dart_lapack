      SUBROUTINE DSCAL(N,DA,DX,INCX);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           DA;
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double           DX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,M,MP1,NINCX;
      // .. Parameters ..
      double           ONE;
      const     ONE=1.0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      if (N <= 0 || INCX <= 0 || DA == ONE) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1


         // clean-up loop

         M = MOD(N,5);
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               DX(I) = DA*DX(I);
            }
            if (N < 5) RETURN;
         }
         MP1 = M + 1;
         DO I = MP1,N,5;
            DX(I) = DA*DX(I);
            DX(I+1) = DA*DX(I+1);
            DX(I+2) = DA*DX(I+2);
            DX(I+3) = DA*DX(I+3);
            DX(I+4) = DA*DX(I+4);
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         DO I = 1,NINCX,INCX;
            DX(I) = DA*DX(I);
         }
      }
      return;

      // End of DSCAL

      }
