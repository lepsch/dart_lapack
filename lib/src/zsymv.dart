// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zsymv(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
  final Complex BETA,
  final Array<Complex> Y_,
  final int INCY,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final Y = Y_.having();
  int I, INFO, IX, IY, J, JX, JY, KX, KY;
  Complex TEMP1, TEMP2;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (LDA < max(1, N)) {
    INFO = 5;
  } else if (INCX == 0) {
    INFO = 7;
  } else if (INCY == 0) {
    INFO = 10;
  }
  if (INFO != 0) {
    xerbla('ZSYMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || ((ALPHA == Complex.zero) && (BETA == Complex.one))) return;

  // Set up the start points in  X  and  Y.

  if (INCX > 0) {
    KX = 1;
  } else {
    KX = 1 - (N - 1) * INCX;
  }
  if (INCY > 0) {
    KY = 1;
  } else {
    KY = 1 - (N - 1) * INCY;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  // First form  y := beta*y.

  if (BETA != Complex.one) {
    if (INCY == 1) {
      if (BETA == Complex.zero) {
        for (I = 1; I <= N; I++) {
          Y[I] = Complex.zero;
        }
      } else {
        for (I = 1; I <= N; I++) {
          Y[I] = BETA * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA == Complex.zero) {
        for (I = 1; I <= N; I++) {
          Y[IY] = Complex.zero;
          IY += INCY;
        }
      } else {
        for (I = 1; I <= N; I++) {
          Y[IY] = BETA * Y[IY];
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == Complex.zero) return;
  if (lsame(UPLO, 'U')) {
    // Form  y  when A is stored in upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        for (I = 1; I <= J - 1; I++) {
          Y[I] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[I];
        }
        Y[J] += TEMP1 * A[J][J] + ALPHA * TEMP2;
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        IX = KX;
        IY = KY;
        for (I = 1; I <= J - 1; I++) {
          Y[IY] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[IX];
          IX += INCX;
          IY += INCY;
        }
        Y[JY] += TEMP1 * A[J][J] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
      }
    }
  } else {
    // Form  y  when A is stored in lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        Y[J] += TEMP1 * A[J][J];
        for (I = J + 1; I <= N; I++) {
          Y[I] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[I];
        }
        Y[J] += ALPHA * TEMP2;
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        Y[JY] += TEMP1 * A[J][J];
        IX = JX;
        IY = JY;
        for (I = J + 1; I <= N; I++) {
          IX += INCX;
          IY += INCY;
          Y[IY] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[IX];
        }
        Y[JY] += ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
      }
    }
  }
}
