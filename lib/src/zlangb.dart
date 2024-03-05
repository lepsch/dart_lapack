import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlangb(
  final String NORM,
  final int N,
  final int KL,
  final int KU,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final WORK = WORK_.having();
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
      // 20
      for (I = max(KU + 2 - J, 1); I <= min(N + KU + 1 - J, KL + KU + 1); I++) {
        // 10
        TEMP = (AB[I][J]).abs();
        if (VALUE < TEMP || disnan(TEMP)) VALUE = TEMP;
      } // 10
    } // 20
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      // 40
      SUM.value = ZERO;
      for (I = max(KU + 2 - J, 1); I <= min(N + KU + 1 - J, KL + KU + 1); I++) {
        // 30
        SUM.value = SUM.value + (AB[I][J]).abs();
      } // 30
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    } // 40
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    for (I = 1; I <= N; I++) {
      // 50
      WORK[I] = ZERO;
    } // 50
    for (J = 1; J <= N; J++) {
      // 70
      K = KU + 1 - J;
      for (I = max(1, J - KU); I <= min(N, J + KL); I++) {
        // 60
        WORK[I] = WORK[I] + (AB[K + I][J]).abs();
      } // 60
    } // 70
    VALUE = ZERO;
    for (I = 1; I <= N; I++) {
      // 80
      TEMP = WORK[I];
      if (VALUE < TEMP || disnan(TEMP)) VALUE = TEMP;
    } // 80
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    for (J = 1; J <= N; J++) {
      // 90
      L = max(1, J - KU);
      K = KU + 1 - J + L;
      zlassq(min(N, J + KL) - L + 1, AB(K, J).asArray(), 1, SCALE, SUM);
    } // 90
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
