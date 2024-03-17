import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/matrix.dart';

double dlangb(
  final String NORM,
  final int N,
  final int KL,
  final int KU,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> WORK,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  const ONE = 1.0, ZERO = 0.0;
  int I, J, K, L;
  double VALUE = 0, TEMP;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      for (I = max(KU + 2 - J, 1); I <= min(N + KU + 1 - J, KL + KU + 1); I++) {
        TEMP = (AB[I][J]).abs();
        if (VALUE < TEMP || disnan(TEMP)) VALUE = TEMP;
      }
    }
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      SUM.value = ZERO;
      for (I = max(KU + 2 - J, 1); I <= min(N + KU + 1 - J, KL + KU + 1); I++) {
        SUM.value += (AB[I][J]).abs();
      }
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    for (I = 1; I <= N; I++) {
      WORK[I] = ZERO;
    }
    for (J = 1; J <= N; J++) {
      K = KU + 1 - J;
      for (I = max(1, J - KU); I <= min(N, J + KL); I++) {
        WORK[I] += (AB[K + I][J]).abs();
      }
    }
    VALUE = ZERO;
    for (I = 1; I <= N; I++) {
      TEMP = WORK[I];
      if (VALUE < TEMP || disnan(TEMP)) VALUE = TEMP;
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    for (J = 1; J <= N; J++) {
      L = max(1, J - KU);
      K = KU + 1 - J + L;
      dlassq(min(N, J + KL) - L + 1, AB(K, J).asArray(), 1, SCALE, SUM);
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
