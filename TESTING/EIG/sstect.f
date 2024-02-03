      SUBROUTINE SSTECT( N, A, B, SHIFT, NUM );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NUM;
      REAL               SHIFT;
      // ..
      // .. Array Arguments ..
      REAL               A( * ), B( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, THREE;
      const              ZERO = 0.0, ONE = 1.0, THREE = 3.0 ;
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

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = SLAMCH( 'Overflow' );

      // Find largest entry

      MX = ABS( A( 1 ) );
      for (I = 1; I <= N - 1; I++) { // 10
         MX = max( MX, ABS( A( I+1 ) ), ABS( B( I ) ) );
      } // 10

      // Handle easy cases, including zero matrix

      if ( SHIFT >= THREE*MX ) {
         NUM = N;
         return;
      }
      if ( SHIFT < -THREE*MX ) {
         NUM = 0;
         return;
      }

      // Compute scale factors as in Kahan's report
      // At this point, MX != 0 so we can divide by it

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

      NUM = 0;
      SSHIFT = ( SHIFT*M1 )*M2;
      U = ( A( 1 )*M1 )*M2 - SSHIFT;
      if ( U <= SUN ) {
         if ( U <= ZERO ) {
            NUM = NUM + 1;
            if (U > -SUN) U = -SUN;
         } else {
            U = SUN;
         }
      }
      for (I = 2; I <= N; I++) { // 20
         TMP = ( B( I-1 )*M1 )*M2;
         U = ( ( A( I )*M1 )*M2-TMP*( TMP / U ) ) - SSHIFT;
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

      // End of SSTECT

      }
