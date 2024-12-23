// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/matrix.dart';

void dlargv(
  final int N,
  final Array<double> X_,
  final int INCX,
  final Array<double> Y_,
  final int INCY,
  final Array<double> C_,
  final int INCC,
) {
  final X = X_.having();
  final Y = Y_.having();
  final C = C_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, IC, IX, IY;
  double F, G, T, TT;

  IX = 1;
  IY = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    F = X[IX];
    G = Y[IY];
    if (G == ZERO) {
      C[IC] = ONE;
    } else if (F == ZERO) {
      C[IC] = ZERO;
      Y[IY] = ONE;
      X[IX] = G;
    } else if (F.abs() > G.abs()) {
      T = G / F;
      TT = sqrt(ONE + T * T);
      C[IC] = ONE / TT;
      Y[IY] = T * C[IC];
      X[IX] = F * TT;
    } else {
      T = F / G;
      TT = sqrt(ONE + T * T);
      Y[IY] = ONE / TT;
      C[IC] = T * Y[IY];
      X[IX] = G * TT;
    }
    IC += INCC;
    IY += INCY;
    IX += INCX;
  }
}
