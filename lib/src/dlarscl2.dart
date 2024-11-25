// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlarscl2(
  final int M,
  final int N,
  final Array<double> D_,
  final Matrix<double> X_,
  final int LDX,
) {
  final D = D_.having();
  final X = X_.having(ld: LDX);
  int I, J;

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      X[I][J] /= D[I];
    }
  }
}
