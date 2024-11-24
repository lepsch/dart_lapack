// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'common.dart';

void dqrt02(
  final int M,
  final int N,
  final int K,
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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

  final EPS = dlamch('Epsilon');

  // Copy the first k columns of the factorization to the array Q

  dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  dlacpy('Lower', M - 1, K, AF(2, 1), LDA, Q(2, 1), LDA);

  // Generate the first n columns of the matrix Q

  srnamc.SRNAMT = 'DORGQR';
  dorgqr(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy R(1:n,1:k)

  dlaset('Full', N, K, ZERO, ZERO, R, LDA);
  dlacpy('Upper', N, K, AF, LDA, R, LDA);

  // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

  dgemm(
      'Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A, LDA, ONE, R, LDA);

  // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = dlange('1', M, K, A, LDA, RWORK);
  var RESID = dlange('1', N, K, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, M)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  dlaset('Full', N, N, ZERO, ONE, R, LDA);
  dsyrk('Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = dlansy('1', 'Upper', N, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, M)) / EPS;
}
