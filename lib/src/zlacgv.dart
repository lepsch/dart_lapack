// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlacgv(
  final int N,
  final Array<Complex> X_,
  final int INCX,
) {
  final X = X_.having();
  if (INCX == 1) {
    for (var I = 1; I <= N; I++) {
      X[I] = X[I].conjugate();
    }
  } else {
    var IOFF = 1;
    if (INCX < 0) IOFF = 1 - (N - 1) * INCX;
    for (var I = 1; I <= N; I++) {
      X[IOFF] = X[IOFF].conjugate();
      IOFF += INCX;
    }
  }
}
