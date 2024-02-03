      SUBROUTINE DLARTGS( X, Y, SIGMA, CS, SN )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      double                  CS, SIGMA, SN, X, Y;
      // ..
*
*  ===================================================================
*
      // .. Parameters ..
      double                  NEGONE, ONE, ZERO;
      PARAMETER          ( NEGONE = -1.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      double                  R, S, THRESH, W, Z;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARTGP
      // ..
      // .. External Functions ..
      double                  DLAMCH;
      // EXTERNAL DLAMCH
      // .. Executable Statements ..
*
      THRESH = DLAMCH('E')
*
      // Compute the first column of B**T*B - SIGMA^2*I, up to a scale
      // factor.
*
      IF( (SIGMA .EQ. ZERO .AND. ABS(X) .LT. THRESH) .OR. (ABS(X) .EQ. SIGMA .AND. Y .EQ. ZERO) ) THEN
         Z = ZERO
         W = ZERO
      ELSE IF( SIGMA .EQ. ZERO ) THEN
         IF( X .GE. ZERO ) THEN
            Z = X
            W = Y
         ELSE
            Z = -X
            W = -Y
         END IF
      ELSE IF( ABS(X) .LT. THRESH ) THEN
         Z = -SIGMA*SIGMA
         W = ZERO
      ELSE
         IF( X .GE. ZERO ) THEN
            S = ONE
         ELSE
            S = NEGONE
         END IF
         Z = S * (ABS(X)-SIGMA) * (S+SIGMA/X)
         W = S * Y
      END IF
*
      // Generate the rotation.
      // CALL DLARTGP( Z, W, CS, SN, R ) might seem more natural;
      // reordering the arguments ensures that if Z = 0 then the rotation
      // is by PI/2.
*
      CALL DLARTGP( W, Z, SN, CS, R )
*
      RETURN
*
      // End DLARTGS
*
      END
