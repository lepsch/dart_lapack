// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgerqf.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorgrq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'common.dart';

void drqt01(
  final int M,
  final int N,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final Matrix<double> Q_,
  final Matrix<double> R_,
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
  final R = R_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e+10;
  final INFO = Box(0);

  final MINMN = min(M, N);
  final EPS = dlamch('Epsilon');

  // Copy the matrix A to the array AF.

  dlacpy('Full', M, N, A, LDA, AF, LDA);

  // Factorize the matrix A in the array AF.

  srnamc.SRNAMT = 'DGERQF';
  dgerqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  if (M <= N) {
    if (M > 0 && M < N) dlacpy('Full', M, N - M, AF, LDA, Q(N - M + 1, 1), LDA);
    if (M > 1) {
      dlacpy('Lower', M - 1, M - 1, AF(2, N - M + 1), LDA,
          Q(N - M + 2, N - M + 1), LDA);
    }
  } else {
    if (N > 1) {
      dlacpy('Lower', N - 1, N - 1, AF(M - N + 2, 1), LDA, Q(2, 1), LDA);
    }
  }

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'DORGRQ';
  dorgrq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy R

  dlaset('Full', M, N, ZERO, ZERO, R, LDA);
  if (M <= N) {
    if (M > 0) {
      dlacpy('Upper', M, M, AF(1, N - M + 1), LDA, R(1, N - M + 1), LDA);
    }
  } else {
    if (M > N && N > 0) dlacpy('Full', M - N, N, AF, LDA, R, LDA);
    if (N > 0) {
      dlacpy('Upper', N, N, AF(M - N + 1, 1), LDA, R(M - N + 1, 1), LDA);
    }
  }

  // Compute R - A*Q'

  dgemm(
      'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA);

  // Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .

  final ANORM = dlange('1', M, N, A, LDA, RWORK);
  var RESID = dlange('1', M, N, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  dlaset('Full', N, N, ZERO, ONE, R, LDA);
  dsyrk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = dlansy('1', 'Upper', N, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
