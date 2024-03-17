import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlansp(
  final String NORM,
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, K;
  double ABSA, VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      K = 1;
      for (J = 1; J <= N; J++) {
        for (I = K; I <= K + J - 1; I++) {
          SUM.value = AP[I].abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        }
        K += J;
      }
    } else {
      K = 1;
      for (J = 1; J <= N; J++) {
        for (I = K; I <= K + N - J; I++) {
          SUM.value = AP[I].abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        }
        K += N - J + 1;
      }
    }
  } else if ((lsame(NORM, 'I')) || (lsame(NORM, 'O')) || (NORM == '1')) {
    // Find normI(A) ( = norm1(A), since A is symmetric).

    VALUE = ZERO;
    K = 1;
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        SUM.value = ZERO;
        for (I = 1; I <= J - 1; I++) {
          ABSA = AP[K].abs();
          SUM.value += ABSA;
          WORK[I] += ABSA;
          K++;
        }
        WORK[J] = SUM.value + AP[K].abs();
        K++;
      }
      for (I = 1; I <= N; I++) {
        SUM.value = WORK[I];
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    } else {
      for (I = 1; I <= N; I++) {
        WORK[I] = ZERO;
      }
      for (J = 1; J <= N; J++) {
        SUM.value = WORK[J] + AP[K].abs();
        K++;
        for (I = J + 1; I <= N; I++) {
          ABSA = AP[K].abs();
          SUM.value += ABSA;
          WORK[I] += ABSA;
          K++;
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    K = 2;
    if (lsame(UPLO, 'U')) {
      for (J = 2; J <= N; J++) {
        zlassq(J - 1, AP(K), 1, SCALE, SUM);
        K += J;
      }
    } else {
      for (J = 1; J <= N - 1; J++) {
        zlassq(N - J, AP(K), 1, SCALE, SUM);
        K += N - J + 1;
      }
    }
    SUM.value = 2 * SUM.value;
    K = 1;
    for (I = 1; I <= N; I++) {
      if (AP[K].toDouble() != ZERO) {
        ABSA = AP[K].toDouble().abs();
        if (SCALE.value < ABSA) {
          SUM.value = ONE + SUM.value * pow(SCALE.value / ABSA, 2);
          SCALE.value = ABSA;
        } else {
          SUM.value += pow(ABSA / SCALE.value, 2);
        }
      }
      if (AP[K].imaginary != ZERO) {
        ABSA = AP[K].imaginary.abs();
        if (SCALE.value < ABSA) {
          SUM.value = ONE + SUM.value * pow(SCALE.value / ABSA, 2);
          SCALE.value = ABSA;
        } else {
          SUM.value += pow(ABSA / SCALE.value, 2);
        }
      }
      if (lsame(UPLO, 'U')) {
        K += I + 1;
      } else {
        K += N - I + 1;
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }
  return VALUE;
}
