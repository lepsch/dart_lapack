      SUBROUTINE ZRSCL( N, A, X, INCX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      COMPLEX*16         A
      // ..
      // .. Array Arguments ..
      COMPLEX*16         X( * )
      // ..

* =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             SAFMAX, SAFMIN, OV, AR, AI, ABSR, ABSI, UR, UI;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      COMPLEX*16         ZLADIV
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

      IF( N.LE.0 ) RETURN

      // Get machine parameters

      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      OV   = DLAMCH( 'O' )

      // Initialize constants related to A.

      AR = DBLE( A )
      AI = DIMAG( A )
      ABSR = ABS( AR )
      ABSI = ABS( AI )

      if ( AI.EQ.ZERO ) {
         // If alpha is real, then we can use csrscl
         zdrscl(N, AR, X, INCX );

      } else if ( AR.EQ.ZERO ) {
         // If alpha has a zero real part, then we follow the same rules as if
         // alpha were real.
         if ( ABSI.GT.SAFMAX ) {
            zdscal(N, SAFMIN, X, INCX );
            zscal(N, DCMPLX( ZERO, -SAFMAX / AI ), X, INCX );
         } else if ( ABSI.LT.SAFMIN ) {
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
         UR = AR + AI * ( AI / AR )
         UI = AI + AR * ( AR / AI )

         if ( (ABS( UR ).LT.SAFMIN).OR.(ABS( UI ).LT.SAFMIN) ) {
            // This means that both alphaR and alphaI are very small.
            zscal(N, DCMPLX( SAFMIN / UR, -SAFMIN / UI ), X, INCX );
            zdscal(N, SAFMAX, X, INCX );
         } else if ( (ABS( UR ).GT.SAFMAX).OR.(ABS( UI ).GT.SAFMAX) ) {
            if ( (ABSR.GT.OV).OR.(ABSI.GT.OV) ) {
               // This means that a and b are both Inf. No need for scaling.
               zscal(N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX );
            } else {
               zdscal(N, SAFMIN, X, INCX );
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
                  zscal(N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX );
               } else {
                  zscal(N, DCMPLX( SAFMAX / UR, -SAFMAX / UI ), X, INCX );
               }
            }
         } else {
            zscal(N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX );
         }
      }

      RETURN

      // End of ZRSCL

      }
