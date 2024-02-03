      SUBROUTINE SLARTGS( X, Y, SIGMA, CS, SN )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL                    CS, SIGMA, SN, X, Y
      // ..

*  ===================================================================

      // .. Parameters ..
      REAL                    NEGONE, ONE, ZERO
      const              NEGONE = -1.0E0, ONE = 1.0E0, ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      REAL                    R, S, THRESH, W, Z
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARTGP
      // ..
      // .. External Functions ..
      REAL                    SLAMCH
      // EXTERNAL SLAMCH
      // .. Executable Statements ..

      THRESH = SLAMCH('E')

      // Compute the first column of B**T*B - SIGMA^2*I, up to a scale
      // factor.

      if ( (SIGMA .EQ. ZERO .AND. ABS(X) .LT. THRESH) .OR. (ABS(X) .EQ. SIGMA .AND. Y .EQ. ZERO) ) {
         Z = ZERO
         W = ZERO
      } else if ( SIGMA .EQ. ZERO ) {
         if ( X .GE. ZERO ) {
            Z = X
            W = Y
         } else {
            Z = -X
            W = -Y
         }
      } else if ( ABS(X) .LT. THRESH ) {
         Z = -SIGMA*SIGMA
         W = ZERO
      } else {
         if ( X .GE. ZERO ) {
            S = ONE
         } else {
            S = NEGONE
         }
         Z = S * (ABS(X)-SIGMA) * (S+SIGMA/X)
         W = S * Y
      }

      // Generate the rotation.
      // CALL SLARTGP( Z, W, CS, SN, R ) might seem more natural;
      // reordering the arguments ensures that if Z = 0 then the rotation
      // is by PI/2.

      slartgp(W, Z, SN, CS, R );

      RETURN

      // End SLARTGS

      }
