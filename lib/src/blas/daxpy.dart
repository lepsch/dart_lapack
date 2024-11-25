// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/blas.dart';

void daxpy(
  final int N,
  final double DA,
  final Array<double> DX_,
  final int INCX,
  final Array<double> DY_,
  final int INCY,
) {
  final DX = DX_.having();
  final DY = DY_.having();

  if (N <= 0) return;
  if (DA == 0.0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    // clean-up loop

    final M = N % 4;
    if (M != 0) {
      for (var I = 1; I <= M; I++) {
        DY[I] += DA * DX[I];
      }
    }
    if (N < 4) return;
    final MP1 = M + 1;
    for (var I = MP1; I <= N; I += 4) {
      DY[I] += DA * DX[I];
      DY[I + 1] += DA * DX[I + 1];
      DY[I + 2] += DA * DX[I + 2];
      DY[I + 3] += DA * DX[I + 3];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    var IX = INCX < 0 ? (-N + 1) * INCX + 1 : 1;
    var IY = INCY < 0 ? (-N + 1) * INCY + 1 : 1;
    for (var I = 1; I <= N; I++) {
      DY[IY] += DA * DX[IX];
      IX += INCX;
      IY += INCY;
    }
  }
}
