// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dger(
  final int M,
  final int N,
  final double ALPHA,
  final Array<double> X_,
  final int INCX,
  final Array<double> Y_,
  final int INCY,
  final Matrix<double> A_,
  final int LDA,
) {
  final X = X_.having();
  final Y = Y_.having();
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
  if (M < 0) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (INCY == 0) {
    INFO = 7;
  } else if (LDA < max(1, M)) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('DGER', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || (ALPHA == ZERO)) return;

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  var JY = INCY > 0 ? 1 : 1 - (N - 1) * INCY;

  if (INCX == 1) {
    for (var J = 1; J <= N; J++) {
      if (Y[JY] != ZERO) {
        final TEMP = ALPHA * Y[JY];
        for (var I = 1; I <= M; I++) {
          A[I][J] += X[I] * TEMP;
        }
      }
      JY += INCY;
    }
  } else {
    final KX = INCX > 0 ? 1 : 1 - (M - 1) * INCX;

    for (var J = 1; J <= N; J++) {
      if (Y[JY] != ZERO) {
        final TEMP = ALPHA * Y[JY];
        var IX = KX;
        for (var I = 1; I <= M; I++) {
          A[I][J] += X[IX] * TEMP;
          IX += INCX;
        }
      }
      JY += INCY;
    }
  }
}
