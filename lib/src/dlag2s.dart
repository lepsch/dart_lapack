// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/slamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlag2s(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> SA_,
  final int LDSA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final SA = SA_.having(ld: LDSA);

  final RMAX = slamch('O');
  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= M; I++) {
      if ((A[I][J] < -RMAX) || (A[I][J] > RMAX)) {
        INFO.value = 1;
        return;
      }
      SA[I][J] = A[I][J];
    }
  }
  INFO.value = 0;
}
