// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dspr(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> X_,
  final int INCX,
  final Array<double> AP,
) {
  final X = X_.having();
  const ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  }
  if (INFO != 0) {
    xerbla('DSPR', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == ZERO)) return;

  // Set the start point in X if the increment is not unity.

  final KX = switch (INCX) {
    <= 0 => 1 - (N - 1) * INCX,
    1 => 0,
    _ => 1,
  };

  // Start the operations. In this version the elements of the array AP
  // are accessed sequentially with one pass through AP.

  var KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  A  when upper triangle is stored in AP.

    if (INCX == 1) {
      for (var J = 1; J <= N; J++) {
        if (X[J] != ZERO) {
          final TEMP = ALPHA * X[J];
          for (var I = 1, K = KK; I <= J; I++, K++) {
            AP[K] += X[I] * TEMP;
          }
        }
        KK += J;
      }
    } else {
      var JX = KX;
      for (var J = 1; J <= N; J++) {
        if (X[JX] != ZERO) {
          final TEMP = ALPHA * X[JX];
          var IX = KX;
          for (var K = KK; K <= KK + J - 1; K++) {
            AP[K] += X[IX] * TEMP;
            IX += INCX;
          }
        }
        JX += INCX;
        KK += J;
      }
    }
  } else {
    // Form  A  when lower triangle is stored in AP.

    if (INCX == 1) {
      for (var J = 1; J <= N; J++) {
        if (X[J] != ZERO) {
          final TEMP = ALPHA * X[J];
          for (var I = J, K = KK; I <= N; I++, K++) {
            AP[K] += X[I] * TEMP;
          }
        }
        KK += N - J + 1;
      }
    } else {
      var JX = KX;
      for (var J = 1; J <= N; J++) {
        if (X[JX] != ZERO) {
          final TEMP = ALPHA * X[JX];
          var IX = JX;
          for (var K = KK; K <= KK + N - J; K++) {
            AP[K] += X[IX] * TEMP;
            IX += INCX;
          }
        }
        JX += INCX;
        KK += N - J + 1;
      }
    }
  }
}
