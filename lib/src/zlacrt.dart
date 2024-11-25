// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlacrt(
  final int N,
  final Array<Complex> CX_,
  final int INCX,
  final Array<Complex> CY_,
  final int INCY,
  final Complex C,
  final Complex S,
) {
  final CX = CX_.having();
  final CY = CY_.having();
  int I, IX, IY;
  Complex CTEMP;

  if (N <= 0) return;
  if (INCX != 1 || INCY != 1) {
    // Code for unequal increments or equal increments not equal to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      CTEMP = C * CX[IX] + S * CY[IY];
      CY[IY] = C * CY[IY] - S * CX[IX];
      CX[IX] = CTEMP;
      IX += INCX;
      IY += INCY;
    }
    return;
  }

  // Code for both increments equal to 1

  for (I = 1; I <= N; I++) {
    CTEMP = C * CX[I] + S * CY[I];
    CY[I] = C * CY[I] - S * CX[I];
    CX[I] = CTEMP;
  }
}
