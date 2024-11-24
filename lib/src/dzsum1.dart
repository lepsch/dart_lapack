// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

double dzsum1(final int N, final Array<Complex> CX_, final int INCX) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final CX = CX_.having();
  int I, NINCX;
  double STEMP;

  STEMP = 0.0;
  if (N <= 0) return 0;
  if (INCX != 1) {
    // CODE FOR INCREMENT NOT EQUAL TO 1
    NINCX = N * INCX;
    for (I = 1; I <= NINCX; I += INCX) {
      STEMP += CX[I].abs();
    }
    return STEMP;
  }

  // CODE FOR INCREMENT EQUAL TO 1
  for (I = 1; I <= N; I++) {
    STEMP += CX[I].abs();
  }
  return STEMP;
}
