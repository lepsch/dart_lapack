      double dlapy2(X, Y ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             X, Y;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             W, XABS, YABS, Z, HUGEVAL;
      bool               X_IS_NAN, Y_IS_NAN;
      // ..
      // .. External Functions ..
      bool               DISNAN;
      // EXTERNAL DISNAN
      // ..
      // .. External Subroutines ..
      double             DLAMCH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      X_IS_NAN = DISNAN( X );
      Y_IS_NAN = DISNAN( Y );
      if (X_IS_NAN) DLAPY2 = X;
      if (Y_IS_NAN) DLAPY2 = Y;
      HUGEVAL = DLAMCH( 'Overflow' );

      if ( !( X_IS_NAN || Y_IS_NAN ) ) {
         XABS = ( X ).abs();
         YABS = ( Y ).abs();
         W = max( XABS, YABS );
         Z = min( XABS, YABS );
         if ( Z == ZERO || W > HUGEVAL ) {
            DLAPY2 = W;
         } else {
            DLAPY2 = W*sqrt( ONE+( Z / W )**2 );
         }
      }
      return;
      }
