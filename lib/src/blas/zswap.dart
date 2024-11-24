// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zswap(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
  final Array<Complex> ZY_,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  final ZY = ZY_.having();
  Complex ZTEMP;
  int I, IX, IY;

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1
    for (I = 1; I <= N; I++) {
      ZTEMP = ZX[I];
      ZX[I] = ZY[I];
      ZY[I] = ZTEMP;
    }
  } else {
    // code for unequal increments or equal increments not equal
    // to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      ZTEMP = ZX[IX];
      ZX[IX] = ZY[IY];
      ZY[IY] = ZTEMP;
      IX += INCX;
      IY += INCY;
    }
  }
}
