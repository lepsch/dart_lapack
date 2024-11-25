// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void drot(
  final int N,
  final Array<double> DX_,
  final int INCX,
  final Array<double> DY_,
  final int INCY,
  final double C,
  final double S,
) {
  final DX = DX_.having();
  final DY = DY_.having();

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (var I = 1; I <= N; I++) {
      final DTEMP = C * DX[I] + S * DY[I];
      DY[I] = C * DY[I] - S * DX[I];
      DX[I] = DTEMP;
    }
  } else {
    // code for unequal increments or equal increments not equal
    // to 1

    var IX = INCX < 0 ? (-N + 1) * INCX + 1 : 1;
    var IY = INCY < 0 ? (-N + 1) * INCY + 1 : 1;
    for (var I = 1; I <= N; I++) {
      final DTEMP = C * DX[IX] + S * DY[IY];
      DY[IY] = C * DY[IY] - S * DX[IX];
      DX[IX] = DTEMP;
      IX += INCX;
      IY += INCY;
    }
  }
}
