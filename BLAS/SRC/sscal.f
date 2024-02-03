      SUBROUTINE SSCAL(N,SA,SX,INCX)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL SA
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int     I,M,MP1,NINCX;
      // ..
      // .. Parameters ..
      REAL ONE
      const     ONE=1.0E+0;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      IF (N.LE.0 .OR. INCX.LE.0 .OR. SA.EQ.ONE) RETURN
      IF (INCX.EQ.1) THEN

         // code for increment equal to 1


         // clean-up loop

         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I) = SA*SX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            SX(I) = SA*SX(I)
            SX(I+1) = SA*SX(I+1)
            SX(I+2) = SA*SX(I+2)
            SX(I+3) = SA*SX(I+3)
            SX(I+4) = SA*SX(I+4)
         END DO
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            SX(I) = SA*SX(I)
         END DO
      END IF
      RETURN

      // End of SSCAL

      }
