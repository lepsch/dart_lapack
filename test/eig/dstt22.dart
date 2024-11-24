// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dstt22(
  final int N,
  final int M,
  final int KBAND,
  final Array<double> AD_,
  final Array<double> AE_,
  final Array<double> SD_,
  final Array<double> SE_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> WORK_,
  final int LDWORK,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AD = AD_.having();
  final AE = AE_.having();
  final SD = SD_.having();
  final SE = SE_.having();
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having(ld: LDWORK);
  final RESULT = RESULT_.having();
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
    ANORM = AD[1].abs() + AE[1].abs();
    for (J = 2; J <= N - 1; J++) {
      ANORM = max(ANORM, AD[J].abs() + AE[J].abs() + AE[J - 1].abs());
    }
    ANORM = max(ANORM, AD[N].abs() + AE[N - 1].abs());
  } else {
    ANORM = AD[1].abs();
  }
  ANORM = max(ANORM, UNFL);

  // Norm of U'AU - S

  for (I = 1; I <= M; I++) {
    for (J = 1; J <= M; J++) {
      WORK[I][J] = ZERO;
      for (K = 1; K <= N; K++) {
        AUKJ = AD[K] * U[K][J];
        if (K != N) AUKJ += AE[K] * U[K + 1][J];
        if (K != 1) AUKJ += AE[K - 1] * U[K - 1][J];
        WORK[I][J] += U[K][I] * AUKJ;
      }
    }
    WORK[I][I] -= SD[I];
    if (KBAND == 1) {
      if (I != 1) WORK[I][I - 1] -= SE[I - 1];
      if (I != N) WORK[I][I + 1] -= SE[I];
    }
  }

  WNORM = dlansy('1', 'L', M, WORK, M, WORK(1, M + 1).asArray());

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (M * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, M * ANORM) / ANORM) / (M * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, M) / (M * ULP);
    }
  }

  // Do Test 2

  // Compute  U'U - I

  dgemm('T', 'N', M, M, N, ONE, U, LDU, U, LDU, ZERO, WORK, M);

  for (J = 1; J <= M; J++) {
    WORK[J][J] -= ONE;
  }

  RESULT[2] =
      min(M, dlange('1', M, M, WORK, M, WORK(1, M + 1).asArray())) / (M * ULP);
}
