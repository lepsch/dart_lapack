      REAL             FUNCTION SGET06( RCOND, RCONDC )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               RCOND, RCONDC
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               EPS, RAT
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' )
      if ( RCOND > ZERO ) {
         if ( RCONDC > ZERO ) {
            RAT = MAX( RCOND, RCONDC ) / MIN( RCOND, RCONDC ) - ( ONE-EPS )
         } else {
            RAT = RCOND / EPS
         }
      } else {
         if ( RCONDC > ZERO ) {
            RAT = RCONDC / EPS
         } else {
            RAT = ZERO
         }
      }
      SGET06 = RAT
      RETURN

      // End of SGET06

      }
