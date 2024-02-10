      double slapy3(final int X, final int Y, final int Z) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               X, Y, Z;
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      double               W, XABS, YABS, ZABS, HUGEVAL;
      // ..
      // .. External Subroutines ..
      double               SLAMCH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT

      HUGEVAL = SLAMCH( 'Overflow' );
      XABS = ( X ).abs();
      YABS = ( Y ).abs();
      ZABS = ( Z ).abs();
      W = max( XABS, YABS, ZABS );
      if ( W == ZERO || W > HUGEVAL ) {
      // W can be zero for max(0,nan,0)
      // adding all three entries together will make sure
      // NaN will not disappear.
         SLAPY3 =  XABS + YABS + ZABS;
      } else {
         SLAPY3 = W*sqrt( ( XABS / W )**2+( YABS / W )**2+ ( ZABS / W )**2 );
      }
      }
