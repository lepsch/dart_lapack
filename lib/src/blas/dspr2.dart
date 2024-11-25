// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dspr2(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> X_,
  final int INCX,
  final Array<double> Y_,
  final int INCY,
  final Array<double> AP_,
) {
  final X = X_.having();
  final Y = Y_.having();
  final AP = AP_.having();
  const ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (INCY == 0) {
    INFO = 7;
  }
  if (INFO != 0) {
    xerbla('DSPR2', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == ZERO)) return;

  // Set up the start points in X and Y if the increments are not both
  // unity.

  int JX = 0, JY = 0, KX = 0, KY = 0;
  if ((INCX != 1) || (INCY != 1)) {
    KX = INCX > 0 ? 1 : 1 - (N - 1) * INCX;
    KY = INCY > 0 ? 1 : 1 - (N - 1) * INCY;
    JX = KX;
    JY = KY;
  }

  // Start the operations. In this version the elements of the array AP
  // are accessed sequentially with one pass through AP.

  var KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  A  when upper triangle is stored in AP.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        if ((X[J] != ZERO) || (Y[J] != ZERO)) {
          final TEMP1 = ALPHA * Y[J];
          final TEMP2 = ALPHA * X[J];
          for (var I = 1, K = KK; I <= J; I++, K++) {
            AP[K] += X[I] * TEMP1 + Y[I] * TEMP2;
          }
        }
        KK += J;
      }
    } else {
      for (var J = 1; J <= N; J++) {
        if ((X[JX] != ZERO) || (Y[JY] != ZERO)) {
          final TEMP1 = ALPHA * Y[JY];
          final TEMP2 = ALPHA * X[JX];
          var IX = KX, IY = KY;
          for (var K = KK; K <= KK + J - 1; K++) {
            AP[K] += X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX += INCX;
            IY += INCY;
          }
        }
        JX += INCX;
        JY += INCY;
        KK += J;
      }
    }
  } else {
    // Form  A  when lower triangle is stored in AP.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        if ((X[J] != ZERO) || (Y[J] != ZERO)) {
          final TEMP1 = ALPHA * Y[J];
          final TEMP2 = ALPHA * X[J];
          for (var I = J, K = KK; I <= N; I++, K++) {
            AP[K] += X[I] * TEMP1 + Y[I] * TEMP2;
          }
        }
        KK += N - J + 1;
      }
    } else {
      for (var J = 1; J <= N; J++) {
        if ((X[JX] != ZERO) || (Y[JY] != ZERO)) {
          final TEMP1 = ALPHA * Y[JY];
          final TEMP2 = ALPHA * X[JX];
          var IX = JX, IY = JY;
          for (var K = KK; K <= KK + N - J; K++) {
            AP[K] += X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX += INCX;
            IY += INCY;
          }
        }
        JX += INCX;
        JY += INCY;
        KK += N - J + 1;
      }
    }
  }
}
