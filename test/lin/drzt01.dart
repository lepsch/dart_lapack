// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dormrz.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

double drzt01(
  final int M,
  final int N,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);

  if (LWORK < M * N + M) {
    xerbla('DRZT01', 8);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  final NORMA = dlange('One-norm', M, N, A, LDA, RWORK);

  // Copy upper triangle R

  dlaset('Full', M, N, ZERO, ZERO, WORK.asMatrix(), M);
  for (var J = 1; J <= M; J++) {
    for (var I = 1; I <= J; I++) {
      WORK[(J - 1) * M + I] = AF[I][J];
    }
  }

  // R *= P(1) * ... *P(m)

  final INFO = Box(0);
  dormrz('Right', 'No transpose', M, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(),
      M, WORK(M * N + 1), LWORK - M * N, INFO);

  // R -= A

  for (var I = 1; I <= N; I++) {
    daxpy(M, -ONE, A(1, I).asArray(), 1, WORK((I - 1) * M + 1), 1);
  }

  var result = dlange('One-norm', M, N, WORK.asMatrix(), M, RWORK);

  result /= dlamch('Epsilon') * max(M, N);
  return NORMA != ZERO ? result / NORMA : result;
}
