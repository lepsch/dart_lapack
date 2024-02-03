      REAL             FUNCTION SCSUM1( N, CX, INCX );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            CX( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, NINCX;
      REAL               STEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      SCSUM1 = 0.0;
      STEMP = 0.0;
      if (N <= 0) RETURN       IF( INCX == 1 ) GO TO 20;

      // CODE FOR INCREMENT NOT EQUAL TO 1

      NINCX = N*INCX;
      DO 10 I = 1, NINCX, INCX;

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ABS( CX( I ) );
      } // 10
      SCSUM1 = STEMP;
      return;

      // CODE FOR INCREMENT EQUAL TO 1

      } // 20
      for (I = 1; I <= N; I++) { // 30

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ABS( CX( I ) );
      } // 30
      SCSUM1 = STEMP;
      return;

      // End of SCSUM1

      }
