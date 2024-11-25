// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/blas.dart';

double dasum(final int N, Array<double> DX, final int INCX) {
  if (N <= 0 || INCX <= 0) return 0.0;

  var DTEMP = 0.0;
  if (INCX == 1) {
    // code for increment equal to 1

    // clean-up loop

    final M = N % 6;
    if (M != 0) {
      for (var I = 1; I <= M; I++) {
        DTEMP += DX[I].abs();
      }
      if (N < 6) {
        return DTEMP;
      }
    }
    final MP1 = M + 1;
    for (var I = MP1; I <= N; I += 6) {
      DTEMP += DX[I].abs() +
          DX[I + 1].abs() +
          DX[I + 2].abs() +
          DX[I + 3].abs() +
          DX[I + 4].abs() +
          DX[I + 5].abs();
    }
  } else {
    // code for increment not equal to 1

    final NINCX = N * INCX;
    for (var I = 1; I <= NINCX; I += INCX) {
      DTEMP += DX[I].abs();
    }
  }
  return DTEMP;
}
