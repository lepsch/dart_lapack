// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeqlf.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorgql.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'common.dart';

void dqlt01(
  final int M,
  final int N,
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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

  final MINMN = min(M, N);
  final EPS = dlamch('Epsilon');

  // Copy the matrix A to the array AF.

  dlacpy('Full', M, N, A, LDA, AF, LDA);

  // Factorize the matrix A in the array AF.

  srnamc.SRNAMT = 'DGEQLF';
  dgeqlf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  dlaset('Full', M, M, ROGUE, ROGUE, Q, LDA);
  if (M >= N) {
    if (N < M && N > 0) dlacpy('Full', M - N, N, AF, LDA, Q(1, M - N + 1), LDA);
    if (N > 1) {
      dlacpy('Upper', N - 1, N - 1, AF(M - N + 1, 2), LDA,
          Q(M - N + 1, M - N + 2), LDA);
    }
  } else {
    if (M > 1) {
      dlacpy('Upper', M - 1, M - 1, AF(1, N - M + 2), LDA, Q(1, 2), LDA);
    }
  }

  // Generate the m-by-m matrix Q

  srnamc.SRNAMT = 'DORGQL';
  dorgql(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy L

  dlaset('Full', M, N, ZERO, ZERO, L, LDA);
  if (M >= N) {
    if (N > 0) {
      dlacpy('Lower', N, N, AF(M - N + 1, 1), LDA, L(M - N + 1, 1), LDA);
    }
  } else {
    if (N > M && M > 0) dlacpy('Full', M, N - M, AF, LDA, L, LDA);
    if (M > 0) {
      dlacpy('Lower', M, M, AF(1, N - M + 1), LDA, L(1, N - M + 1), LDA);
    }
  }

  // Compute L - Q'*A

  dgemm(
      'Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A, LDA, ONE, L, LDA);

  // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = dlange('1', M, N, A, LDA, RWORK);
  var RESID = dlange('1', M, N, L, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, M)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  dlaset('Full', M, M, ZERO, ONE, L, LDA);
  dsyrk('Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = dlansy('1', 'Upper', M, L, LDA, RWORK);

  RESULT[2] = (RESID / max(1, M)) / EPS;
}
