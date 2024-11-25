// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';

import 'zunt01.dart';

void zhst01(
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final int LDH,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDH);
  final Q = Q_.having(ld: LDQ);
  final RESULT = RESULT_.having(length: 2);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0;
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
  zlacpy(' ', N, N, A, LDA, WORK.asMatrix(), LDWORK);

  // Compute Q*H

  zgemm('No transpose', 'No transpose', N, N, N, Complex.one, Q, LDQ, H, LDH,
      Complex.zero, WORK(LDWORK * N + 1).asMatrix(), LDWORK);

  // Compute A - Q*H*Q'

  zgemm(
      'No transpose',
      'Conjugate transpose',
      N,
      N,
      N,
      -Complex.one,
      WORK(LDWORK * N + 1).asMatrix(),
      LDWORK,
      Q,
      LDQ,
      Complex.one,
      WORK.asMatrix(),
      LDWORK);

  ANORM = max(zlange('1', N, N, A, LDA, RWORK), UNFL);
  WNORM = zlange('1', N, N, WORK.asMatrix(), LDWORK, RWORK);

  // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

  RESULT[1] = min(WNORM, ANORM) / max(SMLNUM, ANORM * EPS) / N;

  // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

  zunt01('Columns', N, N, Q, LDQ, WORK, LWORK, RWORK, RESULT(2));
}
