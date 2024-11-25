// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/matrix.dart';

double dla_gerpvgrw(
  final int N,
  final int NCOLS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  int I, J;
  double AMAX, UMAX, RPVGRW;

  RPVGRW = 1.0;

  for (J = 1; J <= NCOLS; J++) {
    AMAX = 0.0;
    UMAX = 0.0;
    for (I = 1; I <= N; I++) {
      AMAX = max(A[I][J].abs(), AMAX);
    }
    for (I = 1; I <= J; I++) {
      UMAX = max(AF[I][J].abs(), UMAX);
    }
    if (UMAX != 0.0) {
      RPVGRW = min(AMAX / UMAX, RPVGRW);
    }
  }
  return RPVGRW;
}
