double sget06(RCOND, RCONDC) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double RCOND, RCONDC;
  // ..

// =====================================================================

  // .. Parameters ..
  double ZERO, ONE;
  const ZERO = 0.0, ONE = 1.0;
  // ..
  // .. Local Scalars ..
  double EPS, RAT;
  // ..
  // .. External Functions ..
  //- REAL               SLAMCH;
  // EXTERNAL SLAMCH
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC MAX, MIN
  // ..
  // .. Executable Statements ..

  EPS = SLAMCH('Epsilon');
  if (RCOND > ZERO) {
    if (RCONDC > ZERO) {
      RAT = max(RCOND, RCONDC) / min(RCOND, RCONDC) - (ONE - EPS);
    } else {
      RAT = RCOND / EPS;
    }
  } else {
    if (RCONDC > ZERO) {
      RAT = RCONDC / EPS;
    } else {
      RAT = ZERO;
    }
  }
  SGET06 = RAT;
  return;
}
