// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/ilatrans.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dla_gbamv(
  final int TRANS,
  final int M,
  final int N,
  final int KL,
  final int KU,
  final double ALPHA,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> X_,
  final int INCX,
  final double BETA,
  final Array<double> Y_,
  final int INCY,
) {
  final AB = AB_.having(ld: LDAB);
  final X = X_.having();
  final Y = Y_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool SYMB_ZERO;
  double TEMP, SAFE1;
  int I, INFO, IY, J, JX, KX, KY, LENX, LENY, KD, KE;

  // Test the input parameters.

  INFO = 0;
  if (!((TRANS == ilatrans('N')) ||
      (TRANS == ilatrans('T')) ||
      (TRANS == ilatrans('C')))) {
    INFO = 1;
  } else if (M < 0) {
    INFO = 2;
  } else if (N < 0) {
    INFO = 3;
  } else if (KL < 0 || KL > M - 1) {
    INFO = 4;
  } else if (KU < 0 || KU > N - 1) {
    INFO = 5;
  } else if (LDAB < KL + KU + 1) {
    INFO = 6;
  } else if (INCX == 0) {
    INFO = 8;
  } else if (INCY == 0) {
    INFO = 11;
  }
  if (INFO != 0) {
    xerbla('DLA_GBAMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  // up the start points in  X  and  Y.

  if (TRANS == ilatrans('N')) {
    LENX = N;
    LENY = M;
  } else {
    LENX = M;
    LENY = N;
  }
  if (INCX > 0) {
    KX = 1;
  } else {
    KX = 1 - (LENX - 1) * INCX;
  }
  if (INCY > 0) {
    KY = 1;
  } else {
    KY = 1 - (LENY - 1) * INCY;
  }

  // Set SAFE1 essentially to be the underflow threshold times the
  // number of additions in each row.

  SAFE1 = dlamch('Safe minimum');
  SAFE1 = (N + 1) * SAFE1;

  // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

  // The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
  // the inexact flag.  Still doesn't help change the iteration order
  // to per-column.

  KD = KU + 1;
  KE = KL + 1;
  IY = KY;
  if (INCX == 1) {
    if (TRANS == ilatrans('N')) {
      for (I = 1; I <= LENY; I++) {
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
          for (J = max(I - KL, 1); J <= min(I + KU, LENX); J++) {
            TEMP = AB[KD + I - J][J].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[J].abs() * TEMP;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);
        IY += INCY;
      }
    } else {
      for (I = 1; I <= LENY; I++) {
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
          for (J = max(I - KL, 1); J <= min(I + KU, LENX); J++) {
            TEMP = AB[KE - I + J][I].abs();
            SYMB_ZERO = SYMB_ZERO && (X[J] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[J].abs() * TEMP;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);
        IY += INCY;
      }
    }
  } else {
    if (TRANS == ilatrans('N')) {
      for (I = 1; I <= LENY; I++) {
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
          JX = KX;
          for (J = max(I - KL, 1); J <= min(I + KU, LENX); J++) {
            TEMP = AB[KD + I - J][J].abs();
            SYMB_ZERO = SYMB_ZERO && (X[JX] == ZERO || TEMP == ZERO);

            Y[IY] += ALPHA * X[JX].abs() * TEMP;
            JX += INCX;
          }
        }
        if (!SYMB_ZERO) Y[IY] += sign(SAFE1, Y[IY]);

        IY += INCY;
      }
    } else {
      for (I = 1; I <= LENY; I++) {
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
          JX = KX;
          for (J = max(I - KL, 1); J <= min(I + KU, LENX); J++) {
            TEMP = AB[KE - I + J][I].abs();
            SYMB_ZERO = SYMB_ZERO && (X[JX] == ZERO || TEMP == ZERO);

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
