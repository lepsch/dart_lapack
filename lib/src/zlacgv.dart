// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlacgv(
  final int N,
  final Array<Complex> X_,
  final int INCX,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  int I, IOFF;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC DCONJG

  if (INCX == 1) {
    for (I = 1; I <= N; I++) {
      X[I] = X[I].conjugate();
    }
  } else {
    IOFF = 1;
    if (INCX < 0) IOFF = 1 - (N - 1) * INCX;
    for (I = 1; I <= N; I++) {
      X[IOFF] = X[IOFF].conjugate();
      IOFF += INCX;
    }
  }
}
