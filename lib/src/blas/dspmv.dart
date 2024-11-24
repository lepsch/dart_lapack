// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dspmv(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> AP_,
  final Array<double> X_,
  final int INCX,
  final double BETA,
  final Array<double> Y_,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final X = X_.having();
  final Y = Y_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 6;
  } else if (INCY == 0) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('DSPMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // Set up the start points in  X  and  Y.

  final KX = INCX > 0 ? 1 : 1 - (N - 1) * INCX;
  final KY = INCY > 0 ? 1 : 1 - (N - 1) * INCY;

  // Start the operations. In this version the elements of the array AP
  // are accessed sequentially with one pass through AP.

  // First form  y := beta*y.

  if (BETA != ONE) {
    if (INCY == 1) {
      if (BETA == ZERO) {
        for (var I = 1; I <= N; I++) {
          Y[I] = ZERO;
        }
      } else {
        for (var I = 1; I <= N; I++) {
          Y[I] *= BETA;
        }
      }
    } else {
      var IY = KY;
      if (BETA == ZERO) {
        for (var I = 1; I <= N; I++) {
          Y[IY] = ZERO;
          IY += INCY;
        }
      } else {
        for (var I = 1; I <= N; I++) {
          Y[IY] *= BETA;
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == ZERO) return;
  var KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  y  when AP contains the upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[J];
        var TEMP2 = ZERO;
        for (var I = 1, K = KK; I <= J - 1; I++, K++) {
          Y[I] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[I];
        }
        Y[J] += TEMP1 * AP[KK + J - 1] + ALPHA * TEMP2;
        KK += J;
      }
    } else {
      var JX = KX;
      var JY = KY;
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[JX];
        var TEMP2 = ZERO;
        var IX = KX, IY = KY;
        for (var K = KK; K <= KK + J - 2; K++) {
          Y[IY] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[IX];
          IX += INCX;
          IY += INCY;
        }
        Y[JY] += TEMP1 * AP[KK + J - 1] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
        KK += J;
      }
    }
  } else {
    // Form  y  when AP contains the lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[J];
        var TEMP2 = ZERO;
        Y[J] += TEMP1 * AP[KK];
        for (var I = J + 1, K = KK + 1; I <= N; I++, K++) {
          Y[I] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[I];
        }
        Y[J] += ALPHA * TEMP2;
        KK += (N - J + 1);
      }
    } else {
      var JX = KX, JY = KY;
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[JX];
        var TEMP2 = ZERO;
        Y[JY] += TEMP1 * AP[KK];
        var IX = JX, IY = JY;
        for (var K = KK + 1; K <= KK + N - J; K++) {
          IX += INCX;
          IY += INCY;
          Y[IY] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[IX];
        }
        Y[JY] += ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
        KK += N - J + 1;
      }
    }
  }
}
