// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zdrot(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
  final Array<Complex> ZY_,
  final int INCY,
  final double C,
  final double S,
) {
  final ZX = ZX_.having();
  final ZY = ZY_.having();

  if (N <= 0) return;

  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (var I = 1; I <= N; I++) {
      final CTEMP = C.toComplex() * ZX[I] + S.toComplex() * ZY[I];
      ZY[I] = C.toComplex() * ZY[I] - S.toComplex() * ZX[I];
      ZX[I] = CTEMP;
    }
  } else {
    // code for unequal increments or equal increments not equal
    // to 1

    var IX = INCX < 0 ? (-N + 1) * INCX + 1 : 1;
    var IY = INCY < 0 ? (-N + 1) * INCY + 1 : 1;
    for (var I = 1; I <= N; I++) {
      final CTEMP = C.toComplex() * ZX[IX] + S.toComplex() * ZY[IY];
      ZY[IY] = C.toComplex() * ZY[IY] - S.toComplex() * ZX[IX];
      ZX[IX] = CTEMP;
      IX += INCX;
      IY += INCY;
    }
  }
}
