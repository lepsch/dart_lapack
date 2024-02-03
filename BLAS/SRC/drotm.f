      SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,INCY,N;
      // ..
      // .. Array Arguments ..
      double           DPARAM(5),DX(*),DY(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double           DFLAG,DH11,DH12,DH21,DH22,TWO,W,Z,ZERO;
      int     I,KX,KY,NSTEPS;
      // ..
      // .. Data statements ..
      final (ZERO,TWO) = (0.0,2.0);
      // ..

      DFLAG = DPARAM(1);
      if (N <= 0 || (DFLAG+TWO == ZERO)) return;
      if (INCX == INCY && INCX > 0) {

         NSTEPS = N*INCX;
         if (DFLAG < ZERO) {
            DH11 = DPARAM(2);
            DH12 = DPARAM(4);
            DH21 = DPARAM(3);
            DH22 = DPARAM(5);
            DO I = 1,NSTEPS,INCX;
               W = DX(I);
               Z = DY(I);
               DX(I) = W*DH11 + Z*DH12;
               DY(I) = W*DH21 + Z*DH22;
            }
         } else if (DFLAG == ZERO) {
            DH12 = DPARAM(4);
            DH21 = DPARAM(3);
            DO I = 1,NSTEPS,INCX;
               W = DX(I);
               Z = DY(I);
               DX(I) = W + Z*DH12;
               DY(I) = W*DH21 + Z;
            }
         } else {
            DH11 = DPARAM(2);
            DH22 = DPARAM(5);
            DO I = 1,NSTEPS,INCX;
               W = DX(I);
               Z = DY(I);
               DX(I) = W*DH11 + Z;
               DY(I) = -W + DH22*Z;
            }
         }
      } else {
         KX = 1;
         KY = 1;
         if (INCX < 0) KX = 1 + (1-N)*INCX;
         if (INCY < 0) KY = 1 + (1-N)*INCY;

         if (DFLAG < ZERO) {
            DH11 = DPARAM(2);
            DH12 = DPARAM(4);
            DH21 = DPARAM(3);
            DH22 = DPARAM(5);
            for (I = 1; I <= N; I++) {
               W = DX(KX);
               Z = DY(KY);
               DX(KX) = W*DH11 + Z*DH12;
               DY(KY) = W*DH21 + Z*DH22;
               KX = KX + INCX;
               KY = KY + INCY;
            }
         } else if (DFLAG == ZERO) {
            DH12 = DPARAM(4);
            DH21 = DPARAM(3);
            for (I = 1; I <= N; I++) {
               W = DX(KX);
               Z = DY(KY);
               DX(KX) = W + Z*DH12;
               DY(KY) = W*DH21 + Z;
               KX = KX + INCX;
               KY = KY + INCY;
            }
         } else {
             DH11 = DPARAM(2);
             DH22 = DPARAM(5);
             for (I = 1; I <= N; I++) {
                W = DX(KX);
                Z = DY(KY);
                DX(KX) = W*DH11 + Z;
                DY(KY) = -W + DH22*Z;
                KX = KX + INCX;
                KY = KY + INCY;
            }
         }
      }
      return;

      // End of DROTM

      }
