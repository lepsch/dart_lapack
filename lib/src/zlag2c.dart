// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/slamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlag2c(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> SA_,
  final int LDSA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final SA = SA_.having(ld: LDSA);
  int I, J;
  double RMAX;

  RMAX = slamch('O');
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      if ((A[I][J].real < -RMAX) ||
          (A[I][J].real > RMAX) ||
          (A[I][J].imaginary < -RMAX) ||
          (A[I][J].imaginary > RMAX)) {
        INFO.value = 1;
        return;
      }
      SA[I][J] = A[I][J];
    }
  }
  INFO.value = 0;
}
