import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/intrinsics/sign.dart';

bool zlctes(final Complex Z, final Complex D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const ZERO = 0.0, ONE = 1.0;

  if (D == Complex.zero) {
    return Z.real < ZERO;
  }

  if (Z.real == ZERO || D.real == ZERO) {
    return sign(ONE, Z.imaginary) != sign(ONE, D.imaginary);
  }

  if (Z.imaginary == ZERO || D.imaginary == ZERO) {
    return sign(ONE, Z.real) != sign(ONE, D.real);
  }

  final ZMAX = max((Z.real).abs(), Z.imaginary.abs());
  return (Z.real / ZMAX) * D.real + (Z.imaginary / ZMAX) * D.imaginary < ZERO;
}
