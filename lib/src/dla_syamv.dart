// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/ilauplo.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dla_syamv(
  final int UPLO,
  final int N,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> X_,
  final int INCX,
  final double BETA,
  final Array<double> Y_,
  final int INCY,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final Y = Y_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool SYMB_ZERO;
  double TEMP, SAFE1;
  int I, INFO, IY, J, JX, KX, KY;

  // Test the input parameters.

  INFO = 0;
  if (UPLO != ilauplo('U') && UPLO != ilauplo('L')) {
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
    xerbla('DLA_SYAMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

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

  // Set SAFE1 essentially to be the underflow threshold times the
  // number of additions in each row.

  SAFE1 = dlamch('Safe minimum');
  SAFE1 = (N + 1) * SAFE1;

  // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

  // The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
  // the inexact flag.  Still doesn't help change the iteration order
  // to per-column.

  IY = KY;
  if (INCX == 1) {
    if (UPLO == ilauplo('U')) {
      for (I = 1; I <= N; I++) {
        if (BETA == ZERO) {
          SYMB_ZERO = true;
          Y[IY] = 0.0;
        } else if (Y[IY] == ZERO) {
          SYMB_ZERO = true;
        } else {
          SYMB_ZERO = false;
          Y[IY] = BETA * Y[IY].abs();
        }
        if (ALPHA != ZERO) {
          for (J = 1; J <= I; J++) {
            TEMP = A[J][I].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[J].abs() * TEMP;
          }
          for (J = I + 1; J <= N; J++) {
            TEMP = A[I][J].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[J].abs() * TEMP;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);

        IY += INCY;
      }
    } else {
      for (I = 1; I <= N; I++) {
        if (BETA == ZERO) {
          SYMB_ZERO = true;
          Y[IY] = 0.0;
        } else if (Y[IY] == ZERO) {
          SYMB_ZERO = true;
        } else {
          SYMB_ZERO = false;
          Y[IY] = BETA * Y[IY].abs();
        }
        if (ALPHA != ZERO) {
          for (J = 1; J <= I; J++) {
            TEMP = A[I][J].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[J].abs() * TEMP;
          }
          for (J = I + 1; J <= N; J++) {
            TEMP = A[J][I].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[J].abs() * TEMP;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);

        IY += INCY;
      }
    }
  } else {
    if (UPLO == ilauplo('U')) {
      for (I = 1; I <= N; I++) {
        if (BETA == ZERO) {
          SYMB_ZERO = true;
          Y[IY] = 0.0;
        } else if (Y[IY] == ZERO) {
          SYMB_ZERO = true;
        } else {
          SYMB_ZERO = false;
          Y[IY] = BETA * Y[IY].abs();
        }
        JX = KX;
        if (ALPHA != ZERO) {
          for (J = 1; J <= I; J++) {
            TEMP = A[J][I].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[JX].abs() * TEMP;
            JX += INCX;
          }
          for (J = I + 1; J <= N; J++) {
            TEMP = A[I][J].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[JX].abs() * TEMP;
            JX += INCX;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);

        IY += INCY;
      }
    } else {
      for (I = 1; I <= N; I++) {
        if (BETA == ZERO) {
          SYMB_ZERO = true;
          Y[IY] = 0.0;
        } else if (Y[IY] == ZERO) {
          SYMB_ZERO = true;
        } else {
          SYMB_ZERO = false;
          Y[IY] = BETA * Y[IY].abs();
        }
        JX = KX;
        if (ALPHA != ZERO) {
          for (J = 1; J <= I; J++) {
            TEMP = A[I][J].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[JX].abs() * TEMP;
            JX += INCX;
          }
          for (J = I + 1; J <= N; J++) {
            TEMP = A[J][I].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[JX].abs() * TEMP;
            JX += INCX;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);

        IY += INCY;
      }
    }
  }
}
