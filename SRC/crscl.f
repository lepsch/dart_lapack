      SUBROUTINE CRSCL( N, A, X, INCX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      COMPLEX            A
      // ..
      // .. Array Arguments ..
      COMPLEX            X( * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               SAFMAX, SAFMIN, OV, AR, AI, ABSR, ABSI, UR , UI
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      COMPLEX            CLADIV
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

      if (N.LE.0) RETURN;

      // Get machine parameters

      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      OV   = SLAMCH( 'O' )

      // Initialize constants related to A.

      AR = REAL( A )
      AI = AIMAG( A )
      ABSR = ABS( AR )
      ABSI = ABS( AI )

      if ( AI == ZERO ) {
         // If alpha is real, then we can use csrscl
         csrscl(N, AR, X, INCX );

      } else if ( AR == ZERO ) {
         // If alpha has a zero real part, then we follow the same rules as if
         // alpha were real.
         if ( ABSI.GT.SAFMAX ) {
            csscal(N, SAFMIN, X, INCX );
            cscal(N, CMPLX( ZERO, -SAFMAX / AI ), X, INCX );
         } else if ( ABSI.LT.SAFMIN ) {
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
         UR = AR + AI * ( AI / AR )
         UI = AI + AR * ( AR / AI )

         if ( (ABS( UR ).LT.SAFMIN).OR.(ABS( UI ).LT.SAFMIN) ) {
            // This means that both alphaR and alphaI are very small.
            cscal(N, CMPLX( SAFMIN / UR, -SAFMIN / UI ), X, INCX );
            csscal(N, SAFMAX, X, INCX );
         } else if ( (ABS( UR ).GT.SAFMAX).OR.(ABS( UI ).GT.SAFMAX) ) {
            if ( (ABSR.GT.OV).OR.(ABSI.GT.OV) ) {
               // This means that a and b are both Inf. No need for scaling.
               cscal(N, CMPLX( ONE / UR, -ONE / UI ), X, INCX );
            } else {
               csscal(N, SAFMIN, X, INCX );
               if ( (ABS( UR ).GT.OV).OR.(ABS( UI ).GT.OV) ) {
                  // Infs were generated. We do proper scaling to avoid them.
                  if ( ABSR.GE.ABSI ) {
                     // ABS( UR ) <= ABS( UI )
                     UR = (SAFMIN * AR) + SAFMIN * (AI * ( AI / AR ))
                     UI = (SAFMIN * AI) + AR * ( (SAFMIN * AR) / AI )
                  } else {
                     // ABS( UR ) > ABS( UI )
                     UR = (SAFMIN * AR) + AI * ( (SAFMIN * AI) / AR )
                     UI = (SAFMIN * AI) + SAFMIN * (AR * ( AR / AI ))
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

      RETURN

      // End of CRSCL

      }
