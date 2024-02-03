      void zdrscl(N, SA, SX, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      double             SA;
      // ..
      // .. Array Arguments ..
      Complex         SX( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
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
      // EXTERNAL ZDSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N <= 0) return;

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;

      // Initialize the denominator to SA and the numerator to 1.

      CDEN = SA;
      CNUM = ONE;

      } // 10
      CDEN1 = CDEN*SMLNUM;
      CNUM1 = CNUM / BIGNUM;
      if ( ( CDEN1 ).abs() > ( CNUM ).abs() && CNUM != ZERO ) {

         // Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.

         MUL = SMLNUM;
         DONE = false;
         CDEN = CDEN1;
      } else if ( ( CNUM1 ).abs() > ( CDEN ).abs() ) {

         // Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.

         MUL = BIGNUM;
         DONE = false;
         CNUM = CNUM1;
      } else {

         // Multiply X by CNUM / CDEN and return.

         MUL = CNUM / CDEN;
         DONE = true;
      }

      // Scale the vector X by MUL

      zdscal(N, MUL, SX, INCX );

      if ( !DONE) GO TO 10;

      return;
      }
