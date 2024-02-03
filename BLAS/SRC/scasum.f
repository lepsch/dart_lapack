      REAL FUNCTION SCASUM(N,CX,INCX)

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      COMPLEX CX(*)
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL STEMP
      int     I,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS,AIMAG,REAL
      // ..
      SCASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      if (INCX.EQ.1) {

         // code for increment equal to 1

         DO I = 1,N
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
         END DO
      } else {

         // code for increment not equal to 1

         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
         END DO
      }
      SCASUM = STEMP
      RETURN

      // End of SCASUM

      }
