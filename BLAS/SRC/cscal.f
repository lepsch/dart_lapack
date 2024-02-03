      SUBROUTINE CSCAL(N,CA,CX,INCX)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX CA
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int     I,NINCX;
      // ..
      // .. Parameters ..
      COMPLEX ONE
      const     ONE= (1.0E+0,0.0E+0);
      // ..
      if (N.LE.0 || INCX.LE.0 || CA == ONE) RETURN;
      if (INCX == 1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            CX(I) = CA*CX(I)
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            CX(I) = CA*CX(I)
         }
      }
      RETURN

      // End of CSCAL

      }
