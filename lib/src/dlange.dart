import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/matrix.dart';

double dlange(
  final String NORM,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final WORK = WORK_.dim(LDA);
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double VALUE = 0, TEMP;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (min(M, N) == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A[i][j])).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= M; I++) {
        TEMP = A[I][J].abs();
        if (VALUE < TEMP || disnan(TEMP)) VALUE = TEMP;
      }
    }
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      SUM.value = ZERO;
      for (I = 1; I <= M; I++) {
        SUM.value = SUM.value + A[I][J].abs();
      }
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    for (I = 1; I <= M; I++) {
      WORK[I] = ZERO;
    }
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I] = WORK[I] + A[I][J].abs();
      }
    }
    VALUE = ZERO;
    for (I = 1; I <= M; I++) {
      TEMP = WORK[I];
      if (VALUE < TEMP || disnan(TEMP)) VALUE = TEMP;
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    for (J = 1; J <= N; J++) {
      dlassq(M, A(1, J).asArray(), 1, SCALE, SUM);
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
