// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zher(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Matrix<Complex> A_,
  final int LDA,
) {
  const ZERO = 0.0;
  final X = X_.having();
  final A = A_.having(ld: LDA);
  Complex TEMP;
  int I, INFO, IX, J, JX, KX = 0;

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
    xerbla('ZHER', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == ZERO)) return;

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
          TEMP = ALPHA.toComplex() * X[J].conjugate();
          for (I = 1; I <= J - 1; I++) {
            A[I][J] += X[I] * TEMP;
          }
          A[J][J] = A[J][J].real.toComplex() + (X[J] * TEMP).real.toComplex();
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[JX].conjugate();
          IX = KX;
          for (I = 1; I <= J - 1; I++) {
            A[I][J] += X[IX] * TEMP;
            IX += INCX;
          }
          A[J][J] = A[J][J].real.toComplex() + (X[JX] * TEMP).real.toComplex();
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
        JX += INCX;
      }
    }
  } else {
    // Form  A  when A is stored in lower triangle.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        if (X[J] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[J].conjugate();
          A[J][J] = A[J][J].real.toComplex() + (TEMP * X[J]).real.toComplex();
          for (I = J + 1; I <= N; I++) {
            A[I][J] += X[I] * TEMP;
          }
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[JX].conjugate();
          A[J][J] = A[J][J].real.toComplex() + (TEMP * X[JX]).real.toComplex();
          IX = JX;
          for (I = J + 1; I <= N; I++) {
            IX += INCX;
            A[I][J] += X[IX] * TEMP;
          }
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
        JX += INCX;
      }
    }
  }
}
