      SUBROUTINE SSTECT( N, A, B, SHIFT, NUM )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NUM;
      REAL               SHIFT
      // ..
      // .. Array Arguments ..
      REAL               A( * ), B( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, THREE
      const              ZERO = 0.0E0, ONE = 1.0E0, THREE = 3.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, TOM, U, UNFL
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )

      // Find largest entry

      MX = ABS( A( 1 ) )
      DO 10 I = 1, N - 1
         MX = MAX( MX, ABS( A( I+1 ) ), ABS( B( I ) ) )
   10 CONTINUE

      // Handle easy cases, including zero matrix

      if ( SHIFT.GE.THREE*MX ) {
         NUM = N
         RETURN
      }
      if ( SHIFT.LT.-THREE*MX ) {
         NUM = 0
         RETURN
      }

      // Compute scale factors as in Kahan's report
      // At this point, MX .NE. 0 so we can divide by it

      SUN = SQRT( UNFL )
      SSUN = SQRT( SUN )
      SOV = SQRT( OVFL )
      TOM = SSUN*SOV
      if ( MX.LE.ONE ) {
         M1 = ONE / MX
         M2 = TOM
      } else {
         M1 = ONE
         M2 = TOM / MX
      }

      // Begin counting

      NUM = 0
      SSHIFT = ( SHIFT*M1 )*M2
      U = ( A( 1 )*M1 )*M2 - SSHIFT
      if ( U.LE.SUN ) {
         if ( U.LE.ZERO ) {
            NUM = NUM + 1
            IF( U.GT.-SUN ) U = -SUN
         } else {
            U = SUN
         }
      }
      DO 20 I = 2, N
         TMP = ( B( I-1 )*M1 )*M2
         U = ( ( A( I )*M1 )*M2-TMP*( TMP / U ) ) - SSHIFT
         if ( U.LE.SUN ) {
            if ( U.LE.ZERO ) {
               NUM = NUM + 1
               IF( U.GT.-SUN ) U = -SUN
            } else {
               U = SUN
            }
         }
   20 CONTINUE
      RETURN

      // End of SSTECT

      }
