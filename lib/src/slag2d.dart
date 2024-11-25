// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';

void slag2d(
  final int M,
  final int N,
  final Matrix<double> SA_,
  final int LDSA,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final SA = SA_.having(ld: LDSA);
  final A = A_.having(ld: LDA);
  int I, J;

  INFO.value = 0;
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      A[I][J] = SA[I][J];
    }
  }
}
