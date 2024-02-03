      REAL             FUNCTION SLAPY3( X, Y, Z );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               X, Y, Z;
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      REAL               W, XABS, YABS, ZABS, HUGEVAL;
      // ..
      // .. External Subroutines ..
      REAL               SLAMCH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      HUGEVAL = SLAMCH( 'Overflow' );
      XABS = ABS( X );
      YABS = ABS( Y );
      ZABS = ABS( Z );
      W = MAX( XABS, YABS, ZABS );
      if ( W == ZERO || W > HUGEVAL ) {
      // W can be zero for max(0,nan,0)
      // adding all three entries together will make sure
      // NaN will not disappear.
         SLAPY3 =  XABS + YABS + ZABS;
      } else {
         SLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+ ( ZABS / W )**2 );
      }
      RETURN;

      // End of SLAPY3

      }
