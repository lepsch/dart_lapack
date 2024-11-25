// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/matrix.dart';

int iladlr(final int M, final int N, final Matrix<double> A_, final int LDA) {
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;

  // Quick test for the common case where one corner is non-zero.
  if (M == 0) return M;

  if (A[M][1] != ZERO || A[M][N] != ZERO) {
    return M;
  }

  // Scan up each column tracking the last zero row seen.
  var result = 0;
  for (var J = 1; J <= N; J++) {
    var I = M;
    while ((A[max(I, 1)][J] == ZERO) && (I >= 1)) {
      I--;
    }
    result = max(result, I);
  }
  return result;
}
