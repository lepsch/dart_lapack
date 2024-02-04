import 'dart:math';

import 'package:lapack/src/install/dlamch.dart';

double dlapy3(final double X, final double Y, final double Z) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  double W, XABS, YABS, ZABS, HUGEVAL;

  HUGEVAL = dlamch('Overflow');
  XABS = (X).abs();
  YABS = (Y).abs();
  ZABS = (Z).abs();
  W = max(XABS, max(YABS, ZABS));
  if (W == ZERO || W > HUGEVAL) {
    // W can be zero for max(0,nan,0)
    // adding all three entries together will make sure
    // NaN will not disappear.
    return XABS + YABS + ZABS;
  }

  return W * sqrt(pow(XABS / W, 2) + pow(YABS / W, 2) + pow(ZABS / W, 2));
}
