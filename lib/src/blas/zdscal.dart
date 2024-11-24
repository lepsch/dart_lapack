// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zdscal(
  final int N,
  final double DA,
  final Array<Complex> ZX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  const ONE = 1.0;

  if (N <= 0 || INCX <= 0 || DA == ONE) return;

  if (INCX == 1) {
    // code for increment equal to 1

    for (var I = 1; I <= N; I++) {
      ZX[I] = Complex(DA * ZX[I].real, DA * ZX[I].imaginary);
    }
  } else {
    // code for increment not equal to 1

    final NINCX = N * INCX;
    for (var I = 1; I <= NINCX; I += INCX) {
      ZX[I] = Complex(DA * ZX[I].real, DA * ZX[I].imaginary);
    }
  }
}
