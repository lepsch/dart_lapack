import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlanhs(
  final String NORM,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      // 20
      for (I = 1; I <= min(N, J + 1); I++) {
        // 10
        SUM.value = A[I][J].abs();
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 10
    } // 20
  } else if ((lsame(NORM, 'O')) || (NORM == '1')) {
    // Find norm1(A).

    VALUE = ZERO;
    for (J = 1; J <= N; J++) {
      // 40
      SUM.value = ZERO;
      for (I = 1; I <= min(N, J + 1); I++) {
        // 30
        SUM.value = SUM.value + (A[I][J]).abs();
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
      for (I = 1; I <= min(N, J + 1); I++) {
        // 60
        WORK[I] += (A[I][J]).abs();
      } // 60
    } // 70
    VALUE = ZERO;
    for (I = 1; I <= N; I++) {
      // 80
      SUM.value = WORK[I];
      if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
    } // 80
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    for (J = 1; J <= N; J++) {
      // 90
      zlassq(min(N, J + 1), A(1, J).asArray(), 1, SCALE, SUM);
    } // 90
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
