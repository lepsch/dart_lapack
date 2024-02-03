      SUBROUTINE DRSCL( N, SA, SX, INCX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      double             SA;
      // ..
      // .. Array Arguments ..
      double             SX( * );
      // ..

* =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               DONE;
      double             BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N.LE.0) RETURN;

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

      // Initialize the denominator to SA and the numerator to 1.

      CDEN = SA
      CNUM = ONE

      } // 10
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      if ( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) {

         // Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.

         MUL = SMLNUM
         DONE = false;
         CDEN = CDEN1
      } else if ( ABS( CNUM1 ).GT.ABS( CDEN ) ) {

         // Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.

         MUL = BIGNUM
         DONE = false;
         CNUM = CNUM1
      } else {

         // Multiply X by CNUM / CDEN and return.

         MUL = CNUM / CDEN
         DONE = true;
      }

      // Scale the vector X by MUL

      dscal(N, MUL, SX, INCX );

      if (.NOT.DONE) GO TO 10;

      RETURN

      // End of DRSCL

      }
