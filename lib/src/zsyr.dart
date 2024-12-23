// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zsyr(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Matrix<Complex> A_,
  final int LDA,
) {
  final X = X_.having();
  final A = A_.having(ld: LDA);
  int I, INFO, IX, J, JX, KX = 0;
  Complex TEMP;

  // Test the input parameters.

  INFO = 0;
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
    xerbla('ZSYR', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == Complex.zero)) return;

  // Set the start point in X if the increment is not unity.

  if (INCX <= 0) {
    KX = 1 - (N - 1) * INCX;
  } else if (INCX != 1) {
    KX = 1;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  if (lsame(UPLO, 'U')) {
    // Form  A  when A is stored in upper triangle.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        if (X[J] != Complex.zero) {
          TEMP = ALPHA * X[J];
          for (I = 1; I <= J; I++) {
            A[I][J] += X[I] * TEMP;
          }
        }
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA * X[JX];
          IX = KX;
          for (I = 1; I <= J; I++) {
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
      for (J = 1; J <= N; J++) {
        if (X[J] != Complex.zero) {
          TEMP = ALPHA * X[J];
          for (I = J; I <= N; I++) {
            A[I][J] += X[I] * TEMP;
          }
        }
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA * X[JX];
          IX = JX;
          for (I = J; I <= N; I++) {
            A[I][J] += X[IX] * TEMP;
            IX += INCX;
          }
        }
        JX += INCX;
      }
    }
  }
}
