// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dget03(
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AINV_,
  final int LDAINV,
  final Matrix<double> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Box<double> RCOND,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final AINV = AINV_.having(ld: LDAINV);
  final WORK = WORK_.having(ld: LDWORK);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlange('1', N, N, A, LDA, RWORK);
  final AINVNM = dlange('1', N, N, AINV, LDAINV, RWORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Compute I - A * AINV

  dgemm('No transpose', 'No transpose', N, N, N, -ONE, AINV, LDAINV, A, LDA,
      ZERO, WORK, LDWORK);
  for (var I = 1; I <= N; I++) {
    WORK[I][I] = ONE + WORK[I][I];
  }

  // Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = dlange('1', N, N, WORK, LDWORK, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / EPS) / N;
}
