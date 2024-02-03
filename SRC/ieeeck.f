      int              FUNCTION IEEECK( ISPEC, ZERO, ONE );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ISPEC;
      REAL               ONE, ZERO;
      // ..

// =====================================================================

      // .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, NEGZRO, NEWZRO, POSINF;
      // ..
      // .. Executable Statements ..
      IEEECK = 1;

      POSINF = ONE / ZERO;
      if ( POSINF <= ONE ) {
         IEEECK = 0;
         return;
      }

      NEGINF = -ONE / ZERO;
      if ( NEGINF >= ZERO ) {
         IEEECK = 0;
         return;
      }

      NEGZRO = ONE / ( NEGINF+ONE );
      if ( NEGZRO != ZERO ) {
         IEEECK = 0;
         return;
      }

      NEGINF = ONE / NEGZRO;
      if ( NEGINF >= ZERO ) {
         IEEECK = 0;
         return;
      }

      NEWZRO = NEGZRO + ZERO;
      if ( NEWZRO != ZERO ) {
         IEEECK = 0;
         return;
      }

      POSINF = ONE / NEWZRO;
      if ( POSINF <= ONE ) {
         IEEECK = 0;
         return;
      }

      NEGINF = NEGINF*POSINF;
      if ( NEGINF >= ZERO ) {
         IEEECK = 0;
         return;
      }

      POSINF = POSINF*POSINF;
      if ( POSINF <= ONE ) {
         IEEECK = 0;
         return;
      }




      // Return if we were only asked to check infinity arithmetic

      if (ISPEC == 0) return;

      NAN1 = POSINF + NEGINF;

      NAN2 = POSINF / NEGINF;

      NAN3 = POSINF / POSINF;

      NAN4 = POSINF*ZERO;

      NAN5 = NEGINF*NEGZRO;

      NAN6 = NAN5*ZERO;

      if ( NAN1 == NAN1 ) {
         IEEECK = 0;
         return;
      }

      if ( NAN2 == NAN2 ) {
         IEEECK = 0;
         return;
      }

      if ( NAN3 == NAN3 ) {
         IEEECK = 0;
         return;
      }

      if ( NAN4 == NAN4 ) {
         IEEECK = 0;
         return;
      }

      if ( NAN5 == NAN5 ) {
         IEEECK = 0;
         return;
      }

      if ( NAN6 == NAN6 ) {
         IEEECK = 0;
         return;
      }

      return;
      }
