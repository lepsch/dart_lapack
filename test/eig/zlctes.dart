bool zlctes(Z, D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex D, Z;
  // ..

// =====================================================================

  // .. Parameters ..

  double ZERO, ONE;
  const ZERO = 0.0, ONE = 1.0;
  Complex CZERO;
  const CZERO = (0.0, 0.0);
  // ..
  // .. Local Scalars ..
  double ZMAX;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS, DBLE, DIMAG, MAX, SIGN
  // ..
  // .. Executable Statements ..

  if (D == CZERO) {
    ZLCTES = (Z.toDouble() < ZERO);
  } else {
    if ((Z).toDouble() == ZERO || D.toDouble() == ZERO) {
      ZLCTES = (sign(ONE, DIMAG(Z)) != sign(ONE, DIMAG(D)));
    } else if (DIMAG(Z) == ZERO || DIMAG(D) == ZERO) {
      ZLCTES = (sign(ONE, (Z).toDouble()) != sign(ONE, D.toDouble()));
    } else {
      ZMAX = max((Z.toDouble()).abs(), (DIMAG(Z)).abs());
      ZLCTES = (((Z).toDouble() / ZMAX) * D.toDouble() +
              (DIMAG(Z) / ZMAX) * DIMAG(D) <
          ZERO);
    }
  }

  return;
}
