// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtrtrs(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ZERO = 0.0, ONE = 1.0;
  bool NOUNIT;

  // Test the input parameters.

  INFO.value = 0;
  NOUNIT = lsame(DIAG, 'N');
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DTRTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check for singularity.

  if (NOUNIT) {
    for (INFO.value = 1; INFO.value <= N; INFO.value++) {
      if (A[INFO.value][INFO.value] == ZERO) return;
    }
  }
  INFO.value = 0;

  // Solve A * x = b  or  A**T * x = b.

  dtrsm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB);
}
