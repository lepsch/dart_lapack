import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/matrix.dart';

double dlanhs(
  final String NORM,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WORK,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    return ZERO;
 }

  if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    var VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= min(N, J + 1); I++) {
        SUM.value = A[I][J].abs();
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    }
    return VALUE;
  }

  if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    var VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      SUM.value = ZERO;
      for (I = 1; I <= min(N, J + 1); I++) {
        SUM.value = SUM.value + A[I][J].abs();
      }
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    }
    return VALUE;
  }

  if (lsame(NORM, 'I')) {
    // Find normI(A).

    for (I = 1; I <= N; I++) {
      WORK[I] = ZERO;
    }
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= min(N, J + 1); I++) {
        WORK[I] = WORK[I] + A[I][J].abs();
      }
    }
    var VALUE = ZERO;
    for (I = 1; I <= N; I++) {
      SUM.value = WORK[I];
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    }
    return VALUE;
  }

  if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    for (J = 1; J <= N; J++) {
      dlassq(min(N, J + 1), A(1, J).asArray(), 1, SCALE, SUM);
    }
    return SCALE.value * sqrt(SUM.value);
  }

  return ZERO;
}
