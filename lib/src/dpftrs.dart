// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dtfsm.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dpftrs(
  final String TRANSR,
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> A_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(offset: zeroIndexedArrayOffset);
  final B = B_.having(ld: LDB, offset: zeroIndexedMatrixOffset);
  const ONE = 1.0;
  bool LOWER, NORMALTRANSR;

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'T')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DPFTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // start execution: there are two triangular solves

  if (LOWER) {
    dtfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB);
    dtfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB);
  } else {
    dtfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB);
    dtfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB);
  }
}
