// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

bool zgennd(final int M, final int N, final Matrix<Complex> A_, final int LDA) {
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;

  final K = min(M, N);
  for (var I = 1; I <= K; I++) {
    final AII = A[I][I];
    if (AII.real < ZERO || AII.imaginary != ZERO) {
      return false;
    }
  }
  return true;
}
