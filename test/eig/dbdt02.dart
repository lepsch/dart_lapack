// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dbdt02(
  final int M,
  final int N,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> C_,
  final int LDC,
  final Matrix<double> U_,
  final int LDU,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having(length: M);
  const ZERO = 0.0, ONE = 1.0;
  int J;
  double BNORM, EPS, REALMN;

  // Quick return if possible

  RESID.value = ZERO;
  if (M <= 0 || N <= 0) return;
  REALMN = (max(M, N)).toDouble();
  EPS = dlamch('Precision');

  // Compute norm(B - U * C)

  for (J = 1; J <= N; J++) {
    dcopy(M, B(1, J).asArray(), 1, WORK, 1);
    dgemv(
        'No transpose', M, M, -ONE, U, LDU, C(1, J).asArray(), 1, ONE, WORK, 1);
    RESID.value = max(RESID.value, dasum(M, WORK, 1));
  }

  // Compute norm of B.

  BNORM = dlange('1', M, N, B, LDB, WORK);

  if (BNORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (BNORM >= RESID.value) {
      RESID.value = (RESID.value / BNORM) / (REALMN * EPS);
    } else {
      if (BNORM < ONE) {
        RESID.value =
            (min(RESID.value, REALMN * BNORM) / BNORM) / (REALMN * EPS);
      } else {
        RESID.value = min(RESID.value / BNORM, REALMN) / (REALMN * EPS);
      }
    }
  }
}
