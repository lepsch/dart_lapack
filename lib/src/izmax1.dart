// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

int izmax1(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  double DMAX;
  int I, IX;

  if (N < 1 || INCX <= 0) return 0;

  if (N == 1) return 1;

  int result = 1;
  if (INCX == 1) {
    // code for increment equal to 1

    DMAX = ZX[1].abs();
    for (I = 2; I <= N; I++) {
      if (ZX[I].abs() > DMAX) {
        result = I;
        DMAX = ZX[I].abs();
      }
    }
  } else {
    // code for increment not equal to 1

    IX = 1;
    DMAX = ZX[1].abs();
    IX += INCX;
    for (I = 2; I <= N; I++) {
      if (ZX[IX].abs() > DMAX) {
        result = I;
        DMAX = ZX[IX].abs();
      }
      IX += INCX;
    }
  }
  return result;
}
