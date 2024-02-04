      double slapy2(X, Y ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double               X, Y;
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      double               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double               W, XABS, YABS, Z, HUGEVAL;
      bool               X_IS_NAN, Y_IS_NAN;
      // ..
      // .. External Functions ..
      //- bool               SISNAN;
      // EXTERNAL SISNAN
      // ..
      // .. External Subroutines ..
      double               SLAMCH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      X_IS_NAN = SISNAN( X );
      Y_IS_NAN = SISNAN( Y );
      if (X_IS_NAN) SLAPY2 = X;
      if (Y_IS_NAN) SLAPY2 = Y;
      HUGEVAL = SLAMCH( 'Overflow' );

      if ( !( X_IS_NAN || Y_IS_NAN ) ) {
         XABS = ( X ).abs();
         YABS = ( Y ).abs();
         W = max( XABS, YABS );
         Z = min( XABS, YABS );
         if ( Z == ZERO || W > HUGEVAL ) {
            SLAPY2 = W;
         } else {
            SLAPY2 = W*sqrt( ONE+( Z / W )**2 );
         }
      }
      return;
      }
