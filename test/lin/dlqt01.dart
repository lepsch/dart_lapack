// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgelqf.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorglq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'common.dart';

void dlqt01(
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

  srnamc.SRNAMT = 'DGELQF';
  dgelqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  if (N > 1) dlacpy('Upper', M, N - 1, AF(1, 2), LDA, Q(1, 2), LDA);

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'DORGLQ';
  dorglq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy L

  dlaset('Full', M, N, ZERO, ZERO, L, LDA);
  dlacpy('Lower', M, N, AF, LDA, L, LDA);

  // Compute L - A*Q'

  dgemm(
      'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, L, LDA);

  // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

  final ANORM = dlange('1', M, N, A, LDA, RWORK);
  final RESID = dlange('1', M, N, L, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  dlaset('Full', N, N, ZERO, ONE, L, LDA);
  dsyrk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  final RESID2 = dlansy('1', 'Upper', N, L, LDA, RWORK);

  RESULT[2] = (RESID2 / max(1, N)) / EPS;
}
