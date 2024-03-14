import 'dart:math';

import 'package:lapack/src/install/dlamch.dart';

double dget06(final double RCOND, final double RCONDC) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;

  final EPS = dlamch('Epsilon');
  if (RCOND > ZERO) {
    if (RCONDC > ZERO) {
      return max(RCOND, RCONDC) / min(RCOND, RCONDC) - (ONE - EPS);
    }

    return RCOND / EPS;
  }

  if (RCONDC > ZERO) {
    return RCONDC / EPS;
  }

  return ZERO;
}
