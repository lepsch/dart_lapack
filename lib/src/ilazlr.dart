// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

int ilazlr(
  final int M,
  final int N,
  final Matrix<Complex> A,
  final int LDA,
) {
  int I, J;

  // Quick test for the common case where one corner is non-zero.
  if (M == 0) {
    return M;
  } else if (A[M][1] != Complex.zero || A[M][N] != Complex.zero) {
    return M;
  } else {
    // Scan up each column tracking the last zero row seen.
    var result = 0;
    for (J = 1; J <= N; J++) {
      I = M;
      while ((A[max(I, 1)][J] == Complex.zero) && (I >= 1)) {
        I--;
      }
      result = max(result, I);
    }
    return result;
  }
}
