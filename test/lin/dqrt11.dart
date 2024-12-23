// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorm2r.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

double dqrt11(
  final int M,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final RDUMMY = Array<double>(1);
  final INFO = Box(0);

  // Test for sufficient workspace

  if (LWORK < M * M + M) {
    xerbla('DQRT11', 7);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0) return ZERO;

  dlaset('Full', M, M, ZERO, ONE, WORK.asMatrix(), M);

  // Form Q

  dorm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK.asMatrix(), M,
      WORK(M * M + 1), INFO);

  // Form Q'*Q

  dorm2r('Left', 'Transpose', M, M, K, A, LDA, TAU, WORK.asMatrix(), M,
      WORK(M * M + 1), INFO);

  for (var J = 1; J <= M; J++) {
    WORK[(J - 1) * M + J] -= ONE;
  }

  return dlange('One-norm', M, M, WORK.asMatrix(), M, RDUMMY) /
      (M * dlamch('Epsilon'));
}
