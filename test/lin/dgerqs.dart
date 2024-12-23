// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dormrq.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgerqs(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || M > N) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LWORK < 1 || LWORK < NRHS && M > 0 && N > 0) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DGERQS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0 || M == 0) return;

  // Solve R*X = B(n-m+1:n,:)

  dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', M, NRHS, ONE,
      A(1, N - M + 1), LDA, B(N - M + 1, 1), LDB);

  // Set B(1:n-m,:) to zero

  dlaset('Full', N - M, NRHS, ZERO, ZERO, B, LDB);

  // B := Q' * B

  dormrq(
      'Left', 'Transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO);
}
