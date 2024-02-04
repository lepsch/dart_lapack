      REAL scsum1(N, CX, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      Complex            CX( * );
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
      if (N <= 0) return;
      IF( INCX == 1 ) GO TO 20;

      // CODE FOR INCREMENT NOT EQUAL TO 1

      NINCX = N*INCX;
      for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) { // 10

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ( CX( I ) ).abs();
      } // 10
      SCSUM1 = STEMP;
      return;

      // CODE FOR INCREMENT EQUAL TO 1

      } // 20
      for (I = 1; I <= N; I++) { // 30

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ( CX( I ) ).abs();
      } // 30
      SCSUM1 = STEMP;
      return;
      }
