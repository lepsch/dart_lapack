// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'dort01.dart';

void dhst01(
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> H_,
  final int LDH,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDH);
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  int LDWORK;
  double ANORM, EPS, SMLNUM, UNFL, WNORM;

  // Quick return if possible

  if (N <= 0) {
    RESULT[1] = ZERO;
    RESULT[2] = ZERO;
    return;
  }

  UNFL = dlamch('Safe minimum');
  EPS = dlamch('Precision');
  SMLNUM = UNFL * N / EPS;

  // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

  // Copy A to WORK

  LDWORK = max(1, N);
  dlacpy(' ', N, N, A, LDA, WORK.asMatrix(LDWORK), LDWORK);

  // Compute Q*H

  dgemm('No transpose', 'No transpose', N, N, N, ONE, Q, LDQ, H, LDH, ZERO,
      WORK(LDWORK * N + 1).asMatrix(LDWORK), LDWORK);

  // Compute A - Q*H*Q'

  dgemm(
      'No transpose',
      'Transpose',
      N,
      N,
      N,
      -ONE,
      WORK(LDWORK * N + 1).asMatrix(LDWORK),
      LDWORK,
      Q,
      LDQ,
      ONE,
      WORK.asMatrix(LDWORK),
      LDWORK);

  ANORM = max(dlange('1', N, N, A, LDA, WORK(LDWORK * N + 1)), UNFL);
  WNORM =
      dlange('1', N, N, WORK.asMatrix(LDWORK), LDWORK, WORK(LDWORK * N + 1));

  // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

  RESULT[1] = min(WNORM, ANORM) / max(SMLNUM, ANORM * EPS) / N;

  // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

  dort01('Columns', N, N, Q, LDQ, WORK, LWORK, RESULT.box(2));
}
