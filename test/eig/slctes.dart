bool slctes(ZR, ZI, D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double D, ZI, ZR;
  // ..

// =====================================================================

  // .. Parameters ..
  double ZERO, ONE;
  const ZERO = 0.0, ONE = 1.0;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC SIGN
  // ..
  // .. Executable Statements ..

  if (D == ZERO) {
    SLCTES = (ZR < ZERO);
  } else {
    SLCTES = (SIGN(ONE, ZR) != SIGN(ONE, D));
  }

  return;
}
