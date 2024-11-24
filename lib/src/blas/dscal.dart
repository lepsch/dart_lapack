// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dscal(
  final int N,
  final double DA,
  final Array<double> DX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();
  const ONE = 1.0;

  if (N <= 0 || INCX <= 0 || DA == ONE) return;

  if (INCX == 1) {
    // code for increment equal to 1

    // clean-up loop

    final M = N % 5;
    if (M != 0) {
      for (var I = 1; I <= M; I++) {
        DX[I] *= DA;
      }
      if (N < 5) return;
    }
    final MP1 = M + 1;
    for (var I = MP1; I <= N; I += 5) {
      DX[I] *= DA;
      DX[I + 1] *= DA;
      DX[I + 2] *= DA;
      DX[I + 3] *= DA;
      DX[I + 4] *= DA;
    }
  } else {
    // code for increment not equal to 1

    final NINCX = N * INCX;
    for (var I = 1; I <= NINCX; I += INCX) {
      DX[I] *= DA;
    }
  }
}
