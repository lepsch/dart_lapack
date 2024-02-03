      SUBROUTINE ZSCAL(N,ZA,ZX,INCX)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ZA
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 ZX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int     I,NINCX;
      // ..
      // .. Parameters ..
      COMPLEX*16 ONE
      const     ONE= (1.0D+0,0.0D+0);
      // ..
      IF (N.LE.0 .OR. INCX.LE.0 .OR. ZA.EQ.ONE) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1

         for (I = 1; I <= N; I++) {
            ZX(I) = ZA*ZX(I)
         }
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = ZA*ZX(I)
         }
      }
      RETURN

      // End of ZSCAL

      }
