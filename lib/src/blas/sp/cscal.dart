      void cscal(N,CA,CX,INCX) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex CA;
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
      Complex ONE;
      const     ONE= (1.0,0.0);
      // ..
      if (N <= 0 || INCX <= 0 || CA == ONE) return;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            CX[I] = CA*CX(I);
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX;
         for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
            CX[I] = CA*CX(I);
         }
      }
      return;
      }
