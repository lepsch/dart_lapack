// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

Complex zdotu(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
  final Array<Complex> ZY_,
  final int INCY,
) {
  final ZX = ZX_.having();
  final ZY = ZY_.having();

  if (N <= 0) return Complex.zero;

  var ZTEMP = Complex.zero;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (var I = 1; I <= N; I++) {
      ZTEMP += ZX[I] * ZY[I];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    var IX = INCX < 0 ? (-N + 1) * INCX + 1 : 1;
    var IY = INCY < 0 ? (-N + 1) * INCY + 1 : 1;
    for (var I = 1; I <= N; I++) {
      ZTEMP += ZX[IX] * ZY[IY];
      IX += INCX;
      IY += INCY;
    }
  }
  return ZTEMP;
}
