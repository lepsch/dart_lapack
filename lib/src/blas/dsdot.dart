// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

double dsdot(
  final int N,
  final Array<double> SX_,
  final int INCX,
  final Array<double> SY_,
  final int INCY,
) {
  final SX = SX_.having();
  final SY = SY_.having();

  var result = 0.0;
  if (N <= 0) return result;
  if (INCX == INCY && INCX > 0) {
    // Code for equal, positive, non-unit increments.

    final NS = N * INCX;
    for (var I = 1; I <= NS; I += INCX) {
      result += SX[I] * SY[I];
    }
  } else {
    // Code for unequal or nonpositive increments.

    var KX = INCX < 0 ? 1 + (1 - N) * INCX : 1;
    var KY = INCY < 0 ? 1 + (1 - N) * INCY : 1;
    for (var I = 1; I <= N; I++) {
      result += SX[KX] * SY[KY];
      KX += INCX;
      KY += INCY;
    }
  }
  return result;
}
