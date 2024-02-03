      double           FUNCTION DLAPY3( X, Y, Z );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             X, Y, Z;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      double             W, XABS, YABS, ZABS, HUGEVAL;
      // ..
      // .. External Subroutines ..
      double             DLAMCH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      HUGEVAL = DLAMCH( 'Overflow' )
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
      // W can be zero for max(0,nan,0)
      // adding all three entries together will make sure
      // NaN will not disappear.
         DLAPY3 =  XABS + YABS + ZABS
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+ ( ZABS / W )**2 )
      END IF
      RETURN

      // End of DLAPY3

      }
