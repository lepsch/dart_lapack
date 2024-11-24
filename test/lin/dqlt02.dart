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
import 'package:lapack/src/dorgql.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'common.dart';

void dqlt02(
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

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) {
    RESULT[1] = ZERO;
    RESULT[2] = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');

  // Copy the last k columns of the factorization to the array Q

  dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  if (K < M) {
    dlacpy('Full', M - K, K, AF(1, N - K + 1), LDA, Q(1, N - K + 1), LDA);
  }
  if (K > 1) {
    dlacpy('Upper', K - 1, K - 1, AF(M - K + 1, N - K + 2), LDA,
        Q(M - K + 1, N - K + 2), LDA);
  }

  // Generate the last n columns of the matrix Q

  srnamc.SRNAMT = 'DORGQL';
  dorgql(M, N, K, Q, LDA, TAU(N - K + 1), WORK, LWORK, INFO);

  // Copy L(m-n+1:m,n-k+1:n)

  dlaset('Full', N, K, ZERO, ZERO, L(M - N + 1, N - K + 1), LDA);
  dlacpy('Lower', K, K, AF(M - K + 1, N - K + 1), LDA, L(M - K + 1, N - K + 1),
      LDA);

  // Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)

  dgemm('Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A(1, N - K + 1),
      LDA, ONE, L(M - N + 1, N - K + 1), LDA);

  // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = dlange('1', M, K, A(1, N - K + 1), LDA, RWORK);
  var RESID = dlange('1', N, K, L(M - N + 1, N - K + 1), LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, M)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  dlaset('Full', N, N, ZERO, ONE, L, LDA);
  dsyrk('Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = dlansy('1', 'Upper', N, L, LDA, RWORK);

  RESULT[2] = (RESID / max(1, M)) / EPS;
}
