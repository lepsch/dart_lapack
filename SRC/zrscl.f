      SUBROUTINE ZRSCL( N, A, X, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INCX, N;
      COMPLEX*16         A
      // ..
      // .. Array Arguments ..
      COMPLEX*16         X( * )
      // ..
*
* =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
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
*
      // Quick return if possible
*
      IF( N.LE.0 ) RETURN
*
      // Get machine parameters
*
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      OV   = DLAMCH( 'O' )
*
      // Initialize constants related to A.
*
      AR = DBLE( A )
      AI = DIMAG( A )
      ABSR = ABS( AR )
      ABSI = ABS( AI )
*
      IF( AI.EQ.ZERO ) THEN
         // If alpha is real, then we can use csrscl
         CALL ZDRSCL( N, AR, X, INCX )
*
      ELSE IF( AR.EQ.ZERO ) THEN
         // If alpha has a zero real part, then we follow the same rules as if
         // alpha were real.
         IF( ABSI.GT.SAFMAX ) THEN
            CALL ZDSCAL( N, SAFMIN, X, INCX )
            CALL ZSCAL( N, DCMPLX( ZERO, -SAFMAX / AI ), X, INCX )
         ELSE IF( ABSI.LT.SAFMIN ) THEN
            CALL ZSCAL( N, DCMPLX( ZERO, -SAFMIN / AI ), X, INCX )
            CALL ZDSCAL( N, SAFMAX, X, INCX )
         ELSE
            CALL ZSCAL( N, DCMPLX( ZERO, -ONE / AI ), X, INCX )
         END IF
*
      ELSE
         // The following numbers can be computed.
         // They are the inverse of the real and imaginary parts of 1/alpha.
         // Note that a and b are always different from zero.
         // NaNs are only possible if either:
         // 1. alphaR or alphaI is NaN.
         // 2. alphaR and alphaI are both infinite, in which case it makes sense
        t // o propagate a NaN.
         UR = AR + AI * ( AI / AR )
         UI = AI + AR * ( AR / AI )
*
         IF( (ABS( UR ).LT.SAFMIN).OR.(ABS( UI ).LT.SAFMIN) ) THEN
            // This means that both alphaR and alphaI are very small.
            CALL ZSCAL( N, DCMPLX( SAFMIN / UR, -SAFMIN / UI ), X, INCX )
            CALL ZDSCAL( N, SAFMAX, X, INCX )
         ELSE IF( (ABS( UR ).GT.SAFMAX).OR.(ABS( UI ).GT.SAFMAX) ) THEN
            IF( (ABSR.GT.OV).OR.(ABSI.GT.OV) ) THEN
               // This means that a and b are both Inf. No need for scaling.
               CALL ZSCAL( N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX )
            ELSE
               CALL ZDSCAL( N, SAFMIN, X, INCX )
               IF( (ABS( UR ).GT.OV).OR.(ABS( UI ).GT.OV) ) THEN
                  // Infs were generated. We do proper scaling to avoid them.
                  IF( ABSR.GE.ABSI ) THEN
                     // ABS( UR ) <= ABS( UI )
                     UR = (SAFMIN * AR) + SAFMIN * (AI * ( AI / AR ))
                     UI = (SAFMIN * AI) + AR * ( (SAFMIN * AR) / AI )
                  ELSE
                     // ABS( UR ) > ABS( UI )
                     UR = (SAFMIN * AR) + AI * ( (SAFMIN * AI) / AR )
                     UI = (SAFMIN * AI) + SAFMIN * (AR * ( AR / AI ))
                  END IF
                  CALL ZSCAL( N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX )
               ELSE
                  CALL ZSCAL( N, DCMPLX( SAFMAX / UR, -SAFMAX / UI ), X, INCX )
               END IF
            END IF
         ELSE
            CALL ZSCAL( N, DCMPLX( ONE / UR, -ONE / UI ), X, INCX )
         END IF
      END IF
*
      RETURN
*
      // End of ZRSCL
*
      END
