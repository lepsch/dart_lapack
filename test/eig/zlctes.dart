import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/intrinsics/sign.dart';

bool zlctes(final Complex Z, final Complex D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const ZERO = 0.0, ONE = 1.0;
  double ZMAX;

  if (D == Complex.zero) {
    return Z.toDouble() < ZERO;
  }

  if ((Z).toDouble() == ZERO || D.toDouble() == ZERO) {
    return (sign(ONE, Z.imaginary) != sign(ONE, D.imaginary));
  } else if (Z.imaginary == ZERO || D.imaginary == ZERO) {
    return (sign(ONE, (Z).toDouble()) != sign(ONE, D.toDouble()));
  } else {
    ZMAX = max((Z.toDouble()).abs(), Z.imaginary.abs());
    return (((Z).toDouble() / ZMAX) * D.toDouble() +
            (Z.imaginary / ZMAX) * D.imaginary <
        ZERO);
  }
}
