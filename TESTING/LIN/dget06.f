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
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
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

      EPS = DLAMCH( 'Epsilon' )
      IF( RCOND.GT.ZERO ) THEN
         IF( RCONDC.GT.ZERO ) THEN
            RAT = MAX( RCOND, RCONDC ) / MIN( RCOND, RCONDC ) - ( ONE-EPS )
         ELSE
            RAT = RCOND / EPS
         END IF
      ELSE
         IF( RCONDC.GT.ZERO ) THEN
            RAT = RCONDC / EPS
         ELSE
            RAT = ZERO
         END IF
      END IF
      DGET06 = RAT
      RETURN

      // End of DGET06

      END
