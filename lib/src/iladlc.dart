// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

int iladlc(final int M, final int N, final Matrix<double> A_, final int LDA) {
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;
  int I;

  // Quick test for the common case where one corner is non-zero.
  if (N == 0) {
    return N;
  }

  if (A[1][N] != ZERO || A[M][N] != ZERO) {
    return N;
  }

  // Now scan each column from the end, returning with the first non-zero.
  int result;
  for (result = N; result >= 1; result--) {
    for (I = 1; I <= M; I++) {
      if (A[I][result] != ZERO) return result;
    }
  }
  return result;
}
