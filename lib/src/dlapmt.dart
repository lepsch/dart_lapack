// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlapmt(
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

  int I, II, IN, J;
  double TEMP;

  if (N <= 1) return;

  for (I = 1; I <= N; I++) {
    K[I] = -K[I];
  }

  if (FORWRD) {
    // Forward permutation

    for (I = 1; I <= N; I++) {
      if (K[I] <= 0) {
        J = I;
        K[J] = -K[J];
        IN = K[J];

        while (K[IN] <= 0) {
          for (II = 1; II <= M; II++) {
            TEMP = X[II][J];
            X[II][J] = X[II][IN];
            X[II][IN] = TEMP;
          }

          K[IN] = -K[IN];
          J = IN;
          IN = K[IN];
        }
      }
    }
  } else {
    // Backward permutation

    for (I = 1; I <= N; I++) {
      if (K[I] <= 0) {
        K[I] = -K[I];
        J = K[I];
        while (J != I) {
          for (II = 1; II <= M; II++) {
            TEMP = X[II][I];
            X[II][I] = X[II][J];
            X[II][J] = TEMP;
          }

          K[J] = -K[J];
          J = K[J];
        }
      }
    }
  }
}
