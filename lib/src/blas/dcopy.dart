// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/blas.dart';

void dcopy(
  final int N,
  final Array<double> DX_,
  final int INCX,
  final Array<double> DY_,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();
  final DY = DY_.having();

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    // clean-up loop

    final M = N % 7;
    if (M != 0) {
      for (var I = 1; I <= M; I++) {
        DY[I] = DX[I];
      }
      if (N < 7) return;
    }
    final MP1 = M + 1;
    for (var I = MP1; I <= N; I += 7) {
      DY[I] = DX[I];
      DY[I + 1] = DX[I + 1];
      DY[I + 2] = DX[I + 2];
      DY[I + 3] = DX[I + 3];
      DY[I + 4] = DX[I + 4];
      DY[I + 5] = DX[I + 5];
      DY[I + 6] = DX[I + 6];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    var IX = 1;
    var IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (var I = 1; I <= N; I++) {
      DY[IY] = DX[IX];
      IX += INCX;
      IY += INCY;
    }
  }
}
