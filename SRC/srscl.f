      SUBROUTINE SRSCL( N, SA, SX, INCX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      REAL               SA
      // ..
      // .. Array Arguments ..
      REAL               SX( * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               DONE;
      REAL               BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

      // Initialize the denominator to SA and the numerator to 1.

      CDEN = SA
      CNUM = ONE

   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      if ( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) {

         // Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.

         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      } else if ( ABS( CNUM1 ).GT.ABS( CDEN ) ) {

         // Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.

         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      } else {

         // Multiply X by CNUM / CDEN and return.

         MUL = CNUM / CDEN
         DONE = .TRUE.
      }

      // Scale the vector X by MUL

      CALL SSCAL( N, MUL, SX, INCX )

      IF( .NOT.DONE ) GO TO 10

      RETURN

      // End of SRSCL

      }
