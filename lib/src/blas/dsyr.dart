// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dsyr(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> X_,
  final int INCX,
  final Matrix<double> A_,
  final int LDA,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (LDA < max(1, N)) {
    INFO = 7;
  }
  if (INFO != 0) {
    xerbla('DSYR', INFO);
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

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  if (lsame(UPLO, 'U')) {
    // Form  A  when A is stored in upper triangle.

    if (INCX == 1) {
      for (var J = 1; J <= N; J++) {
        if (X[J] != ZERO) {
          final TEMP = ALPHA * X[J];
          for (var I = 1; I <= J; I++) {
            A[I][J] += X[I] * TEMP;
          }
        }
      }
    } else {
      var JX = KX;
      for (var J = 1; J <= N; J++) {
        if (X[JX] != ZERO) {
          final TEMP = ALPHA * X[JX];
          var IX = KX;
          for (var I = 1; I <= J; I++) {
            A[I][J] += X[IX] * TEMP;
            IX += INCX;
          }
        }
        JX += INCX;
      }
    }
  } else {
    // Form  A  when A is stored in lower triangle.

    if (INCX == 1) {
      for (var J = 1; J <= N; J++) {
        if (X[J] != ZERO) {
          final TEMP = ALPHA * X[J];
          for (var I = J; I <= N; I++) {
            A[I][J] += X[I] * TEMP;
          }
        }
      }
    } else {
      var JX = KX;
      for (var J = 1; J <= N; J++) {
        if (X[JX] != ZERO) {
          final TEMP = ALPHA * X[JX];
          var IX = JX;
          for (var I = J; I <= N; I++) {
            A[I][J] += X[IX] * TEMP;
            IX += INCX;
          }
        }
        JX += INCX;
      }
    }
  }
}
