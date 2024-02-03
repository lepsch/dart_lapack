      SUBROUTINE DSCAL(N,DA,DX,INCX)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           DA;
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      double           DX(*);
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int     I,M,MP1,NINCX;
      // .. Parameters ..
      double           ONE;
      const     ONE=1.0D+0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      IF (N.LE.0 .OR. INCX.LE.0 .OR. DA.EQ.ONE) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1


         // clean-up loop

         M = MOD(N,5)
         if (M.NE.0) {
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         }
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      }
      RETURN

      // End of DSCAL

      }
