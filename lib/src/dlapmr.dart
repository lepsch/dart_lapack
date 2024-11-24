// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlapmr(
  final bool FORWRD,
  final int M,
  final int N,
  final Matrix<double> X_,
  final int LDX,
  final Array<int> K,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having(ld: LDX);
  int I, IN, J, JJ;
  double TEMP;

  if (M <= 1) return;

  for (I = 1; I <= M; I++) {
    K[I] = -K[I];
  }

  if (FORWRD) {
    // Forward permutation

    for (I = 1; I <= M; I++) {
      if (K[I] > 0) continue;

      J = I;
      K[J] = -K[J];
      IN = K[J];

      while (true) {
        if (K[IN] > 0) break;

        for (JJ = 1; JJ <= N; JJ++) {
          TEMP = X[J][JJ];
          X[J][JJ] = X[IN][JJ];
          X[IN][JJ] = TEMP;
        }

        K[IN] = -K[IN];
        J = IN;
        IN = K[IN];
      }
    }
  } else {
    // Backward permutation

    for (I = 1; I <= M; I++) {
      if (K[I] > 0) continue;

      K[I] = -K[I];
      J = K[I];
      while (true) {
        if (J == I) break;

        for (JJ = 1; JJ <= N; JJ++) {
          TEMP = X[I][JJ];
          X[I][JJ] = X[J][JJ];
          X[J][JJ] = TEMP;
        }

        K[J] = -K[J];
        J = K[J];
      }
    }
  }
}
