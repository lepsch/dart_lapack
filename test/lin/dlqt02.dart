// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorglq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'common.dart';

void dlqt02(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final Matrix<double> Q_,
  final Matrix<double> L_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final L = L_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e+10;
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');

  // Copy the first k rows of the factorization to the array Q

  dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  dlacpy('Upper', K, N - 1, AF(1, 2), LDA, Q(1, 2), LDA);

  // Generate the first n columns of the matrix Q

  srnamc.SRNAMT = 'DORGLQ';
  dorglq(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy L(1:k,1:m)

  dlaset('Full', K, M, ZERO, ZERO, L, LDA);
  dlacpy('Lower', K, M, AF, LDA, L, LDA);

  // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'

  dgemm(
      'No transpose', 'Transpose', K, M, N, -ONE, A, LDA, Q, LDA, ONE, L, LDA);

  // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .

  final ANORM = dlange('1', K, N, A, LDA, RWORK);
  var RESID = dlange('1', K, M, L, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  dlaset('Full', M, M, ZERO, ONE, L, LDA);
  dsyrk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = dlansy('1', 'Upper', M, L, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
