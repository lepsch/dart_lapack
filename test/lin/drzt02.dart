// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dormrz.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

double drzt02(
  final int M,
  final int N,
  final Matrix<double> AF_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
) {
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  if (LWORK < N * N + N) {
    xerbla('DRZT02', 7);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  // Q := I

  dlaset('Full', N, N, ZERO, ONE, WORK.asMatrix(), N);

  // Q := P(1) * ... * P(m) * Q

  dormrz('Left', 'No transpose', N, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(),
      N, WORK(N * N + 1), LWORK - N * N, INFO);

  // Q := P(m) * ... * P(1) * Q

  dormrz('Left', 'Transpose', N, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(), N,
      WORK(N * N + 1), LWORK - N * N, INFO);

  // Q := Q - I

  for (var I = 1; I <= N; I++) {
    WORK[(I - 1) * N + I] -= ONE;
  }

  return dlange('One-norm', N, N, WORK.asMatrix(), N, RWORK) /
      (dlamch('Epsilon') * max(M, N));
}
