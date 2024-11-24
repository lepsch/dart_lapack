// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

int ilazlc(
  final int M,
  final int N,
  final Matrix<Complex> A,
  final int LDA,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // Quick test for the common case where one corner is non-zero.
  if (N == 0) {
    return N;
  } else if (A[1][N] != Complex.zero || A[M][N] != Complex.zero) {
    return N;
  } else {
    // Now scan each column from the end, returning with the first non-zero.
    for (var J = N; J >= 1; J--) {
      for (var I = 1; I <= M; I++) {
        if (A[I][J] != Complex.zero) return J;
      }
    }
  }
  return N;
}
