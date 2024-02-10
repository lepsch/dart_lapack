import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dstt22(
  final int N,
  final int M,
  final int KBAND,
  final Array<double> AD,
  final Array<double> AE,
  final Array<double> SD,
  final Array<double> SE,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> WORK,
  final int LDWORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int I, J, K;
  double ANORM, AUKJ, ULP, UNFL, WNORM;

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0 || M <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon');

  // Do Test 1

  // Compute the 1-norm of A.

  if (N > 1) {
    ANORM = (AD[1]).abs() + (AE[1]).abs();
    for (J = 2; J <= N - 1; J++) {
      ANORM = max(ANORM, (AD[J]).abs() + (AE[J]).abs() + (AE[J - 1]).abs());
    }
    ANORM = max(ANORM, (AD[N]).abs() + (AE[N - 1]).abs());
  } else {
    ANORM = (AD[1]).abs();
  }
  ANORM = max(ANORM, UNFL);

  // Norm of U'AU - S

  for (I = 1; I <= M; I++) {
    for (J = 1; J <= M; J++) {
      WORK[I][J] = ZERO;
      for (K = 1; K <= N; K++) {
        AUKJ = AD[K] * U[K][J];
        if (K != N) AUKJ = AUKJ + AE[K] * U[K + 1][J];
        if (K != 1) AUKJ = AUKJ + AE[K - 1] * U[K - 1][J];
        WORK[I][J] = WORK[I][J] + U[K][I] * AUKJ;
      }
    }
    WORK[I][I] = WORK[I][I] - SD[I];
    if (KBAND == 1) {
      if (I != 1) WORK[I][I - 1] = WORK[I][I - 1] - SE[I - 1];
      if (I != N) WORK[I][I + 1] = WORK[I][I + 1] - SE[I];
    }
  }

  WNORM = dlansy('1', 'L', M, WORK, M, WORK[1][M + 1]);

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (M * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, M * ANORM) / ANORM) / (M * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, M.toDouble()) / (M * ULP);
    }
  }

  // Do Test 2

  // Compute  U'U - I

  dgemm('T', 'N', M, M, N, ONE, U, LDU, U, LDU, ZERO, WORK, M);

  for (J = 1; J <= M; J++) {
    WORK[J][J] = WORK[J][J] - ONE;
  }

  RESULT[2] =
      min(M.toDouble(), dlange('1', M, M, WORK, M, WORK(1, M + 1).asArray())) /
          (M * ULP);
}
