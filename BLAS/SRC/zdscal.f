      SUBROUTINE ZDSCAL(N,DA,ZX,INCX);

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           DA;
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 ZX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,NINCX;
      // .. Parameters ..
      double           ONE;
      const     ONE=1.0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG
      // ..
      if (N <= 0 || INCX <= 0 || DA == ONE) return;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            ZX(I) = DCMPLX(DA*DBLE(ZX(I)),DA*DIMAG(ZX(I)));
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         DO I = 1,NINCX,INCX;
            ZX(I) = DCMPLX(DA*DBLE(ZX(I)),DA*DIMAG(ZX(I)));
         }
      }
      return;

      // End of ZDSCAL

      }
