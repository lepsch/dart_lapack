      void crscl(N, A, X, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      Complex            A;
      // ..
      // .. Array Arguments ..
      Complex            X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double               SAFMAX, SAFMIN, OV, AR, AI, ABSR, ABSI, UR , UI;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- COMPLEX            CLADIV;
      // EXTERNAL SLAMCH, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSSCAL, CSRSCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N <= 0) return;

      // Get machine parameters

      SAFMIN = SLAMCH( 'S' );
      SAFMAX = ONE / SAFMIN;
      OV   = SLAMCH( 'O' );

      // Initialize constants related to A.

      AR = double( A );
      AI = AIMAG( A );
      ABSR = ( AR ).abs();
      ABSI = ( AI ).abs();

      if ( AI == ZERO ) {
         // If alpha is real, then we can use csrscl
         csrscl(N, AR, X, INCX );

      } else if ( AR == ZERO ) {
         // If alpha has a zero real part, then we follow the same rules as if
         // alpha were real.
         if ( ABSI > SAFMAX ) {
            csscal(N, SAFMIN, X, INCX );
            cscal(N, CMPLX( ZERO, -SAFMAX / AI ), X, INCX );
         } else if ( ABSI < SAFMIN ) {
            cscal(N, CMPLX( ZERO, -SAFMIN / AI ), X, INCX );
            csscal(N, SAFMAX, X, INCX );
         } else {
            cscal(N, CMPLX( ZERO, -ONE / AI ), X, INCX );
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
            cscal(N, CMPLX( SAFMIN / UR, -SAFMIN / UI ), X, INCX );
            csscal(N, SAFMAX, X, INCX );
         } else if ( (( UR ).abs() > SAFMAX) || (( UI ).abs() > SAFMAX) ) {
            if ( (ABSR > OV) || (ABSI > OV) ) {
               // This means that a and b are both Inf. No need for scaling.
               cscal(N, CMPLX( ONE / UR, -ONE / UI ), X, INCX );
            } else {
               csscal(N, SAFMIN, X, INCX );
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
                  cscal(N, CMPLX( ONE / UR, -ONE / UI ), X, INCX );
               } else {
                  cscal(N, CMPLX( SAFMAX / UR, -SAFMAX / UI ), X, INCX );
               }
            }
         } else {
            cscal(N, CMPLX( ONE / UR, -ONE / UI ), X, INCX );
         }
      }

      return;
      }