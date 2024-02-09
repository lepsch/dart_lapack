import 'dart:math';

import 'package:lapack/src/matrix.dart';

double dsxt1(
  final int IJOB,
  final Array<double> D1,
  final int N1,
  final Array<double> D2,
  final int N2,
  final double ABSTOL,
  final double ULP,
  final double UNFL,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  int I, J;
  double TEMP1, TEMP2;

  TEMP1 = ZERO;

  J = 1;
  for (I = 1; I <= N1; I++) {
    while (D2[J] < D1[I] && J < N2) {
      J = J + 1;
    }
    if (J == 1) {
      TEMP2 = (D2[J] - D1[I]).abs();
      if (IJOB == 2) TEMP2 = TEMP2 / max(UNFL, ABSTOL + ULP * (D1[I]).abs());
    } else {
      TEMP2 = min((D2[J] - D1[I]).abs(), (D1[I] - D2[J - 1]).abs());
      if (IJOB == 2) TEMP2 = TEMP2 / max(UNFL, ABSTOL + ULP * (D1[I]).abs());
    }
    TEMP1 = max(TEMP1, TEMP2);
  }

  return TEMP1;
}
