// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

int idamax(final int N, final Array<double> DX_, final int INCX) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();

  if (N < 1 || INCX <= 0) return 0;

  if (N == 1) return 1;

  var index = 1;
  if (INCX == 1) {
    // code for increment equal to 1

    var DMAX = DX[1].abs();
    for (var I = 2; I <= N; I++) {
      if (DX[I].abs() > DMAX) {
        index = I;
        DMAX = DX[I].abs();
      }
    }
  } else {
    // code for increment not equal to 1

    var IX = 1;
    var DMAX = DX[1].abs();
    IX += INCX;
    for (var I = 2; I <= N; I++) {
      if (DX[IX].abs() > DMAX) {
        index = I;
        DMAX = DX[IX].abs();
      }
      IX += INCX;
    }
  }
  return index;
}
