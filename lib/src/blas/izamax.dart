// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/blas/dcabs1.dart';
import 'package:dart_lapack/src/matrix.dart';

int izamax(final int N, final Array<Complex> ZX_, final int INCX) {
  final ZX = ZX_.having();

  if (N < 1 || INCX <= 0) return 0;

  if (N == 1) return 1;

  var index = 1;
  if (INCX == 1) {
    // code for increment equal to 1

    var DMAX = dcabs1(ZX[1]);
    for (var I = 2; I <= N; I++) {
      if (dcabs1(ZX[I]) > DMAX) {
        index = I;
        DMAX = dcabs1(ZX[I]);
      }
    }
  } else {
    // code for increment not equal to 1

    var IX = 1;
    var DMAX = dcabs1(ZX[1]);
    IX += INCX;
    for (var I = 2; I <= N; I++) {
      if (dcabs1(ZX[IX]) > DMAX) {
        index = I;
        DMAX = dcabs1(ZX[IX]);
      }
      IX += INCX;
    }
  }
  return index;
}
