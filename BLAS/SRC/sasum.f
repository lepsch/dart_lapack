      REAL FUNCTION SASUM(N,SX,INCX)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL STEMP
      int     I,M,MP1,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS,MOD
      // ..
      SASUM = 0.0;
      STEMP = 0.0;
      if (N <= 0 || INCX <= 0) RETURN;
      if (INCX == 1) {
         // code for increment equal to 1


         // clean-up loop

         M = MOD(N,6)
         if (M != 0) {
            for (I = 1; I <= M; I++) {
               STEMP = STEMP + ABS(SX(I))
            }
            if (N < 6) {
               SASUM = STEMP
               RETURN
            }
         }
         MP1 = M + 1
         DO I = MP1,N,6
            STEMP = STEMP + ABS(SX(I)) + ABS(SX(I+1)) + ABS(SX(I+2)) + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + ABS(SX(I))
         }
      }
      SASUM = STEMP
      RETURN

      // End of SASUM

      }
