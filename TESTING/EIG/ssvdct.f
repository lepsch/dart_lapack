      void ssvdct(N, S, E, SHIFT, NUM ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NUM;
      REAL               SHIFT;
      // ..
      // .. Array Arguments ..
      REAL               E( * ), S( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, TOM, U, UNFL;
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine constants

      UNFL = 2*SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;

      // Find largest entry

      MX = ABS( S( 1 ) );
      for (I = 1; I <= N - 1; I++) { // 10
         MX = max( MX, ABS( S( I+1 ) ), ABS( E( I ) ) );
      } // 10

      if ( MX == ZERO ) {
         if ( SHIFT < ZERO ) {
            NUM = 0;
         } else {
            NUM = 2*N;
         }
         return;
      }

      // Compute scale factors as in Kahan's report

      SUN = sqrt( UNFL );
      SSUN = sqrt( SUN );
      SOV = sqrt( OVFL );
      TOM = SSUN*SOV;
      if ( MX <= ONE ) {
         M1 = ONE / MX;
         M2 = TOM;
      } else {
         M1 = ONE;
         M2 = TOM / MX;
      }

      // Begin counting

      U = ONE;
      NUM = 0;
      SSHIFT = ( SHIFT*M1 )*M2;
      U = -SSHIFT;
      if ( U <= SUN ) {
         if ( U <= ZERO ) {
            NUM = NUM + 1;
            if (U > -SUN) U = -SUN;
         } else {
            U = SUN;
         }
      }
      TMP = ( S( 1 )*M1 )*M2;
      U = -TMP*( TMP / U ) - SSHIFT;
      if ( U <= SUN ) {
         if ( U <= ZERO ) {
            NUM = NUM + 1;
            if (U > -SUN) U = -SUN;
         } else {
            U = SUN;
         }
      }
      for (I = 1; I <= N - 1; I++) { // 20
         TMP = ( E( I )*M1 )*M2;
         U = -TMP*( TMP / U ) - SSHIFT;
         if ( U <= SUN ) {
            if ( U <= ZERO ) {
               NUM = NUM + 1;
               if (U > -SUN) U = -SUN;
            } else {
               U = SUN;
            }
         }
         TMP = ( S( I+1 )*M1 )*M2;
         U = -TMP*( TMP / U ) - SSHIFT;
         if ( U <= SUN ) {
            if ( U <= ZERO ) {
               NUM = NUM + 1;
               if (U > -SUN) U = -SUN;
            } else {
               U = SUN;
            }
         }
      } // 20
      return;
      }
