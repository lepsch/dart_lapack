// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtpttr(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final AP = AP_.having();
  bool LOWER;
  int I, J, K;

  // Test the input parameters.

  INFO.value = 0;
  LOWER = lsame(UPLO, 'L');
  if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DTPTTR', -INFO.value);
    return;
  }

  if (LOWER) {
    K = 0;
    for (J = 1; J <= N; J++) {
      for (I = J; I <= N; I++) {
        K++;
        A[I][J] = AP[K];
      }
    }
  } else {
    K = 0;
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        K++;
        A[I][J] = AP[K];
      }
    }
  }
}
