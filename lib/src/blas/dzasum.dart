// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/blas/dcabs1.dart';
import 'package:dart_lapack/src/matrix.dart';

double dzasum(final int N, final Array<Complex> ZX_, final int INCX) {
  final ZX = ZX_.having();
  if (N <= 0 || INCX <= 0) return 0.0;

  var STEMP = 0.0;
  if (INCX == 1) {
    // code for increment equal to 1

    for (var I = 1; I <= N; I++) {
      STEMP += dcabs1(ZX[I]);
    }
  } else {
    // code for increment not equal to 1

    final NINCX = N * INCX;
    for (var I = 1; I <= NINCX; I += INCX) {
      STEMP += dcabs1(ZX[I]);
    }
  }
  return STEMP;
}
