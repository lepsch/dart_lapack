import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlansy(
  final String NORM,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double ABSA, VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        // 20
        for (I = 1; I <= J; I++) {
          // 10
          SUM.value = (A[I][J]).abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        } // 10
      } // 20
    } else {
      for (J = 1; J <= N; J++) {
        // 40
        for (I = J; I <= N; I++) {
          // 30
          SUM.value = (A[I][J]).abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        } // 30
      } // 40
    }
  } else if ((lsame(NORM, 'I')) || (lsame(NORM, 'O')) || (NORM == '1')) {
    // Find normI(A) ( = norm1(A), since A is symmetric).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        // 60
        SUM.value = ZERO;
        for (I = 1; I <= J - 1; I++) {
          // 50
          ABSA = (A[I][J]).abs();
          SUM.value = SUM.value + ABSA;
          WORK[I] += ABSA;
        } // 50
        WORK[J] = SUM.value + (A[J][J]).abs();
      } // 60
      for (I = 1; I <= N; I++) {
        // 70
        SUM.value = WORK[I];
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 70
    } else {
      for (I = 1; I <= N; I++) {
        // 80
        WORK[I] = ZERO;
      } // 80
      for (J = 1; J <= N; J++) {
        // 100
        SUM.value = WORK[J] + (A[J][J]).abs();
        for (I = J + 1; I <= N; I++) {
          // 90
          ABSA = (A[I][J]).abs();
          SUM.value = SUM.value + ABSA;
          WORK[I] += ABSA;
        } // 90
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 100
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    if (lsame(UPLO, 'U')) {
      for (J = 2; J <= N; J++) {
        // 110
        zlassq(J - 1, A(1, J).asArray(), 1, SCALE, SUM);
      } // 110
    } else {
      for (J = 1; J <= N - 1; J++) {
        // 120
        zlassq(N - J, A(J + 1, J).asArray(), 1, SCALE, SUM);
      } // 120
    }
    SUM.value = 2 * SUM.value;
    zlassq(N, A.asArray(), LDA + 1, SCALE, SUM);
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
