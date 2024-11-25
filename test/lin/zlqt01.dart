// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgelqf.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zunglq.dart';

import 'common.dart';

void zlqt01(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final Matrix<Complex> Q_,
  final Matrix<Complex> L_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
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
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final INFO = Box(0);

  final MINMN = min(M, N);
  final EPS = dlamch('Epsilon');

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N, A, LDA, AF, LDA);

  // Factorize the matrix A in the array AF.

  srnamc.SRNAMT = 'ZGELQF';
  zgelqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  if (N > 1) zlacpy('Upper', M, N - 1, AF(1, 2), LDA, Q(1, 2), LDA);

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'ZUNGLQ';
  zunglq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy L

  zlaset('Full', M, N, Complex.zero, Complex.zero, L, LDA);
  zlacpy('Lower', M, N, AF, LDA, L, LDA);

  // Compute L - A*Q'

  zgemm('No transpose', 'Conjugate transpose', M, N, N, Complex(-ONE), A, LDA,
      Q, LDA, Complex.one, L, LDA);

  // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

  final ANORM = zlange('1', M, N, A, LDA, RWORK);
  var RESID = zlange('1', M, N, L, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  zlaset('Full', N, N, Complex.zero, Complex.one, L, LDA);
  zherk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = zlansy('1', 'Upper', N, L, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
