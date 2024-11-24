// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void icopy(
  final int N,
  final Array<int> SX_,
  final int INCX,
  final Array<int> SY_,
  final int INCY,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SX = SX_.having();
  final SY = SY_.having();
  int MP1;
  if (N <= 0) return;

  if (INCX != 1 || INCY != 1) {
    // Code for unequal increments or equal increments not equal to 1

    var IX = 1;
    var IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (var I = 1; I <= N; I++) {
      SY[IY] = SX[IX];
      IX += INCX;
      IY += INCY;
    }
    return;
  }

  // Code for both increments equal to 1

  // Clean-up loop

  final M = (N % 7);
  if (M != 0) {
    for (var I = 1; I <= M; I++) {
      SY[I] = SX[I];
    }
    if (N < 7) return;
  }
  MP1 = M + 1;
  for (var I = MP1; I <= N; I += 7) {
    SY[I] = SX[I];
    SY[I + 1] = SX[I + 1];
    SY[I + 2] = SX[I + 2];
    SY[I + 3] = SX[I + 3];
    SY[I + 4] = SX[I + 4];
    SY[I + 5] = SX[I + 5];
    SY[I + 6] = SX[I + 6];
  }
}
