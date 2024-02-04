      void zrscl(N, A, X, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      Complex         A;
      // ..
      // .. Array Arguments ..
      Complex         X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             SAFMAX, SAFMIN, OV, AR, AI, ABSR, ABSI, UR, UI;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      //- Complex         ZLADIV;
      // EXTERNAL DLAMCH, ZLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, ZDSCAL, ZDRSCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N <= 0) return;

      // Get machine parameters

      SAFMIN = dlamch( 'S' );
      SAFMAX = ONE / SAFMIN;
      OV   = dlamch( 'O' );

      // Initialize constants related to A.

      AR = A.toDouble();
      AI = DIMAG( A );
      ABSR = ( AR ).abs();
      ABSI = ( AI ).abs();

      if ( AI == ZERO ) {
         // If alpha is real, then we can use csrscl
         zdrscl(N, AR, X, INCX );

      } else if ( AR == ZERO ) {
         // If alpha has a zero real part, then we follow the same rules as if
         // alpha were real.
         if ( ABSI > SAFMAX ) {
            zdscal(N, SAFMIN, X, INCX );
            zscal(N, DCMPLX( ZERO, -SAFMAX / AI ), X, INCX );
         } else if ( ABSI < SAFMIN ) {
            zscal(N, DCMPLX( ZERO, -SAFMIN / AI ), X, INCX );
            zdscal(N, SAFMAX, X, INCX );
         } else {
            zscal(N, DCMPLX( ZERO, -ONE / AI ), X, INCX );
         }

      } else {
         // The following numbers can be computed.
         // They are the inverse of the real and imaginary parts of 1/alpha.
         // Note that a and b are always different from zero.
         // NaNs are only possible if either:
         // 1. alphaR or alphaI is NaN.
         // 2. alphaR and alphaI are both infinite, in which case it makes sense
         // to propagate a NaN.
         UR = AR + AI * ( AI / AR );
         UI = AI + AR * ( AR / AI );

         if ( (( UR ).abs() < SAFMIN) || (( UI ).abs() < SAFMIN) ) {
            // This means that both alphaR and alphaI are very small.
            zscal(N, DCMPLX( SAFMIN / UR, -SAFMIN / UI ), X, INCX );
            zdscal(N, SAFMAX, X, INCX );
         } else if ( (( UR ).abs() > SAFMAX) || (( UI ).abs() > SAFMAX) ) {
            if ( (ABSR > OV) || (ABSI > OV) ) {
               // This means that a and b are both Inf. No need for scaling.
               zscal(N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX );
            } else {
               zdscal(N, SAFMIN, X, INCX );
               if ( (( UR ).abs() > OV) || (( UI ).abs() > OV) ) {
                  // Infs were generated. We do proper scaling to avoid them.
                  if ( ABSR >= ABSI ) {
                     // ABS( UR ) <= ABS( UI )
                     UR = (SAFMIN * AR) + SAFMIN * (AI * ( AI / AR ));
                     UI = (SAFMIN * AI) + AR * ( (SAFMIN * AR) / AI );
                  } else {
                     // ABS( UR ) > ABS( UI )
                     UR = (SAFMIN * AR) + AI * ( (SAFMIN * AI) / AR );
                     UI = (SAFMIN * AI) + SAFMIN * (AR * ( AR / AI ));
                  }
                  zscal(N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX );
               } else {
                  zscal(N, DCMPLX( SAFMAX / UR, -SAFMAX / UI ), X, INCX );
               }
            }
         } else {
            zscal(N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX );
         }
      }

      return;
      }
