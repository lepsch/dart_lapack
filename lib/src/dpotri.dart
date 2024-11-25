// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlauum.dart';
import 'package:dart_lapack/src/dtrtri.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dpotri(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DPOTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  dtrtri(UPLO, 'Non-unit', N, A, LDA, INFO);
  if (INFO.value > 0) return;

  // Form inv(U) * inv(U)**T or inv(L)**T * inv(L).

  dlauum(UPLO, N, A, LDA, INFO);
}
