void slartgs(X, Y, SIGMA, CS, SN) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double CS, SIGMA, SN, X, Y;
  // ..

// ===================================================================

  // .. Parameters ..
  double NEGONE, ONE, ZERO;
  const NEGONE = -1.0, ONE = 1.0, ZERO = 0.0;
  // ..
  // .. Local Scalars ..
  double R, S, THRESH, W, Z;
  // ..
  // .. External Subroutines ..
  // EXTERNAL SLARTGP
  // ..
  // .. External Functions ..
  //- REAL                    SLAMCH;
  // EXTERNAL SLAMCH

  THRESH = SLAMCH('E');

  // Compute the first column of B**T*B - SIGMA^2*I, up to a scale
  // factor.

  if ((SIGMA == ZERO && (X).abs() < THRESH) ||
      ((X).abs() == SIGMA && Y == ZERO)) {
    Z = ZERO;
    W = ZERO;
  } else if (SIGMA == ZERO) {
    if (X >= ZERO) {
      Z = X;
      W = Y;
    } else {
      Z = -X;
      W = -Y;
    }
  } else if ((X).abs() < THRESH) {
    Z = -SIGMA * SIGMA;
    W = ZERO;
  } else {
    if (X >= ZERO) {
      S = ONE;
    } else {
      S = NEGONE;
    }
    Z = S * ((X).abs() - SIGMA) * (S + SIGMA / X);
    W = S * Y;
  }

  // Generate the rotation.
  // CALL SLARTGP( Z, W, CS, SN, R ) might seem more natural;
  // reordering the arguments ensures that if Z = 0 then the rotation
  // is by PI/2.

  slartgp(W, Z, SN, CS, R);

  return;

  // End SLARTGS
}
