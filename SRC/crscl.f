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

      IF( N.LE.0 ) RETURN

      // Get machine parameters

      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      OV   = SLAMCH( 'O' )

      // Initialize constants related to A.

      AR = REAL( A )
      AI = AIMAG( A )
      ABSR = ABS( AR )
      ABSI = ABS( AI )

      IF( AI.EQ.ZERO ) THEN
         // If alpha is real, then we can use csrscl
         CALL CSRSCL( N, AR, X, INCX )

      ELSE IF( AR.EQ.ZERO ) THEN
         // If alpha has a zero real part, then we follow the same rules as if
         // alpha were real.
         IF( ABSI.GT.SAFMAX ) THEN
            CALL CSSCAL( N, SAFMIN, X, INCX )
            CALL CSCAL( N, CMPLX( ZERO, -SAFMAX / AI ), X, INCX )
         ELSE IF( ABSI.LT.SAFMIN ) THEN
            CALL CSCAL( N, CMPLX( ZERO, -SAFMIN / AI ), X, INCX )
            CALL CSSCAL( N, SAFMAX, X, INCX )
         } else {
            CALL CSCAL( N, CMPLX( ZERO, -ONE / AI ), X, INCX )
         END IF

      } else {
         // The following numbers can be computed.
         // They are the inverse of the real and imaginary parts of 1/alpha.
         // Note that a and b are always different from zero.
         // NaNs are only possible if either:
         // 1. alphaR or alphaI is NaN.
         // 2. alphaR and alphaI are both infinite, in which case it makes sense
        t // o propagate a NaN.
         UR = AR + AI * ( AI / AR )
         UI = AI + AR * ( AR / AI )

         IF( (ABS( UR ).LT.SAFMIN).OR.(ABS( UI ).LT.SAFMIN) ) THEN
            // This means that both alphaR and alphaI are very small.
            CALL CSCAL( N, CMPLX( SAFMIN / UR, -SAFMIN / UI ), X, INCX )
            CALL CSSCAL( N, SAFMAX, X, INCX )
         ELSE IF( (ABS( UR ).GT.SAFMAX).OR.(ABS( UI ).GT.SAFMAX) ) THEN
            IF( (ABSR.GT.OV).OR.(ABSI.GT.OV) ) THEN
               // This means that a and b are both Inf. No need for scaling.
               CALL CSCAL( N, CMPLX( ONE / UR, -ONE / UI ), X, INCX )
            } else {
               CALL CSSCAL( N, SAFMIN, X, INCX )
               IF( (ABS( UR ).GT.OV).OR.(ABS( UI ).GT.OV) ) THEN
                  // Infs were generated. We do proper scaling to avoid them.
                  IF( ABSR.GE.ABSI ) THEN
                     // ABS( UR ) <= ABS( UI )
                     UR = (SAFMIN * AR) + SAFMIN * (AI * ( AI / AR ))
                     UI = (SAFMIN * AI) + AR * ( (SAFMIN * AR) / AI )
                  } else {
                     // ABS( UR ) > ABS( UI )
                     UR = (SAFMIN * AR) + AI * ( (SAFMIN * AI) / AR )
                     UI = (SAFMIN * AI) + SAFMIN * (AR * ( AR / AI ))
                  END IF
                  CALL CSCAL( N, CMPLX( ONE / UR, -ONE / UI ), X, INCX )
               } else {
                  CALL CSCAL( N, CMPLX( SAFMAX / UR, -SAFMAX / UI ), X, INCX )
               END IF
            END IF
         } else {
            CALL CSCAL( N, CMPLX( ONE / UR, -ONE / UI ), X, INCX )
         END IF
      END IF

      RETURN

      // End of CRSCL

      }
