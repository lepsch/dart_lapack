      void csscal(N,SA,CX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double SA;
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      Complex CX(*);
      // ..

// =====================================================================

      // .. Local Scalars ..
      int     I,NINCX;
      // ..
      // .. Parameters ..
      double ONE;
      const     ONE=1.0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG,CMPLX,REAL
      // ..
      if (N <= 0 || INCX <= 0 || SA == ONE) return;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            CX[I] = CMPLX(SA*double(CX(I)),SA*AIMAG(CX(I)));
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
            CX[I] = CMPLX(SA*double(CX(I)),SA*AIMAG(CX(I)));
         }
      }
      return;
      }
