import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlanhb(
  final String NORM,
  final String UPLO,
  final int N,
  final int K,
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
  int I, J, L;
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
        for (I = max(K + 2 - J, 1); I <= K; I++) {
          // 10
          SUM.value = (AB[I][J]).abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        } // 10
        SUM.value = AB[K + 1][J].toDouble().abs();
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      } // 20
    } else {
      for (J = 1; J <= N; J++) {
        // 40
        SUM.value = AB[1][J].toDouble().abs();
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        for (I = 2; I <= min(N + 1 - J, K + 1); I++) {
          // 30
          SUM.value = (AB[I][J]).abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        } // 30
      } // 40
    }
  } else if ((lsame(NORM, 'I')) || (lsame(NORM, 'O')) || (NORM == '1')) {
    // Find normI(A) ( = norm1(A), since A is hermitian).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        // 60
        SUM.value = ZERO;
        L = K + 1 - J;
        for (I = max(1, J - K); I <= J - 1; I++) {
          // 50
          ABSA = AB[L + I][J].abs();
          SUM.value = SUM.value + ABSA;
          WORK[I] += ABSA;
        } // 50
        WORK[J] = SUM.value + AB[K + 1][J].toDouble().abs();
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
        SUM.value = WORK[J] + AB[1][J].toDouble().abs();
        L = 1 - J;
        for (I = J + 1; I <= min(N, J + K); I++) {
          // 90
          ABSA = (AB[L + I][J]).abs();
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
    if (K > 0) {
      if (lsame(UPLO, 'U')) {
        for (J = 2; J <= N; J++) {
          // 110
          zlassq(
              min(J - 1, K), AB(max(K + 2 - J, 1), J).asArray(), 1, SCALE, SUM);
        } // 110
        L = K + 1;
      } else {
        for (J = 1; J <= N - 1; J++) {
          // 120
          zlassq(min(N - J, K), AB(2, J).asArray(), 1, SCALE, SUM);
        } // 120
        L = 1;
      }
      SUM.value = 2 * SUM.value;
    } else {
      L = 1;
    }
    for (J = 1; J <= N; J++) {
      // 130
      if (AB[L][J].toDouble() != ZERO) {
        ABSA = AB[L][J].toDouble().abs();
        if (SCALE.value < ABSA) {
          SUM.value = ONE + SUM.value * pow(SCALE.value / ABSA, 2);
          SCALE.value = ABSA;
        } else {
          SUM.value = SUM.value + pow(ABSA / SCALE.value, 2);
        }
      }
    } // 130
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
