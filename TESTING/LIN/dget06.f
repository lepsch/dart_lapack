      double           FUNCTION DGET06( RCOND, RCONDC );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             RCOND, RCONDC;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             EPS, RAT;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' );
      if ( RCOND > ZERO ) {
         if ( RCONDC > ZERO ) {
            RAT = MAX( RCOND, RCONDC ) / MIN( RCOND, RCONDC ) - ( ONE-EPS );
         } else {
            RAT = RCOND / EPS;
         }
      } else {
         if ( RCONDC > ZERO ) {
            RAT = RCONDC / EPS;
         } else {
            RAT = ZERO;
         }
      }
      DGET06 = RAT;
      RETURN;

      // End of DGET06

      }
