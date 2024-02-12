import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void dlarra(
  final int N,
  final Array<double> D,
  final Array<double> E,
  final Array<double> E2,
  final double SPLTOL,
  final double TNRM,
  final Box<int> NSPLIT,
  final Array<int> ISPLIT,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  int I;
  double EABS, TMP1;

  INFO.value = 0;
  NSPLIT.value = 1;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  // Compute splitting points
  if (SPLTOL < ZERO) {
    // Criterion based on absolute off-diagonal value
    TMP1 = (SPLTOL).abs() * TNRM;
    for (I = 1; I <= N - 1; I++) {
      EABS = E[I].abs();
      if (EABS <= TMP1) {
        E[I] = ZERO;
        E2[I] = ZERO;
        ISPLIT[NSPLIT.value] = I;
        NSPLIT.value = NSPLIT.value + 1;
      }
    }
  } else {
    // Criterion that guarantees relative accuracy
    for (I = 1; I <= N - 1; I++) {
      EABS = E[I].abs();
      if (EABS <= SPLTOL * sqrt(D[I].abs()) * sqrt(D[I + 1].abs())) {
        E[I] = ZERO;
        E2[I] = ZERO;
        ISPLIT[NSPLIT.value] = I;
        NSPLIT.value = NSPLIT.value + 1;
      }
    }
  }
  ISPLIT[NSPLIT.value] = N;
}
