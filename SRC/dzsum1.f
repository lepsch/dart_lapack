      double           FUNCTION DZSUM1( N, CX, INCX );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         CX( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, NINCX;
      double             STEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      DZSUM1 = 0.0D0
      STEMP = 0.0D0
      IF( N.LE.0 ) RETURN       IF( INCX.EQ.1 ) GO TO 20

      // CODE FOR INCREMENT NOT EQUAL TO 1

      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      DZSUM1 = STEMP
      RETURN

      // CODE FOR INCREMENT EQUAL TO 1

   20 CONTINUE
      DO 30 I = 1, N

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      DZSUM1 = STEMP
      RETURN

      // End of DZSUM1

      END
