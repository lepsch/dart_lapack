// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zgerc(
  final int M,
  final int N,
  final Complex ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Array<Complex> Y_,
  final int INCY,
  final Matrix<Complex> A_,
  final int LDA,
) {
  final X = X_.having();
  final Y = Y_.having();
  final A = A_.having(ld: LDA);

  Complex TEMP;
  int I, INFO, IX, J, JY, KX;

  // Test the input parameters.

  INFO = 0;
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
    xerbla('ZGERC', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || (ALPHA == Complex.zero)) return;

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  if (INCY > 0) {
    JY = 1;
  } else {
    JY = 1 - (N - 1) * INCY;
  }
  if (INCX == 1) {
    for (J = 1; J <= N; J++) {
      if (Y[JY] != Complex.zero) {
        TEMP = ALPHA * Y[JY].conjugate();
        for (I = 1; I <= M; I++) {
          A[I][J] += X[I] * TEMP;
        }
      }
      JY += INCY;
    }
  } else {
    if (INCX > 0) {
      KX = 1;
    } else {
      KX = 1 - (M - 1) * INCX;
    }
    for (J = 1; J <= N; J++) {
      if (Y[JY] != Complex.zero) {
        TEMP = ALPHA * Y[JY].conjugate();
        IX = KX;
        for (I = 1; I <= M; I++) {
          A[I][J] += X[IX] * TEMP;
          IX += INCX;
        }
      }
      JY += INCY;
    }
  }
}
