// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/ztfsm.dart';

void zpftrs(
  final String TRANSR,
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> A_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(offset: zeroIndexedArrayOffset);
  final B = B_.having(ld: LDB);
  bool LOWER, NORMALTRANSR;

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'C')) {
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
    xerbla('ZPFTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // start execution: there are two triangular solves

  if (LOWER) {
    ztfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, Complex.one, A, B, LDB);
    ztfsm(TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, Complex.one, A, B, LDB);
  } else {
    ztfsm(TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, Complex.one, A, B, LDB);
    ztfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, Complex.one, A, B, LDB);
  }
}
