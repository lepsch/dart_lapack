// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zscal(
  final int N,
  final Complex ZA,
  final Array<Complex> ZX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  int I, NINCX;

  if (N <= 0 || INCX <= 0 || ZA == Complex.one) return;
  if (INCX == 1) {
    // code for increment equal to 1

    for (I = 1; I <= N; I++) {
      ZX[I] = ZA * ZX[I];
    }
  } else {
    // code for increment not equal to 1

    NINCX = N * INCX;
    for (I = 1; I <= NINCX; I += INCX) {
      ZX[I] = ZA * ZX[I];
    }
  }
}
