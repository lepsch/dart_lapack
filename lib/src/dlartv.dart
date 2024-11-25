// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlartv(
  final int N,
  final Array<double> X_,
  final int INCX,
  final Array<double> Y_,
  final int INCY,
  final Array<double> C_,
  final Array<double> S_,
  final int INCC,
) {
  final X = X_.having();
  final Y = Y_.having();
  final C = C_.having();
  final S = S_.having();
  int I, IC, IX, IY;
  double XI, YI;

  IX = 1;
  IY = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    XI = X[IX];
    YI = Y[IY];
    X[IX] = C[IC] * XI + S[IC] * YI;
    Y[IY] = C[IC] * YI - S[IC] * XI;
    IX += INCX;
    IY += INCY;
    IC += INCC;
  }
}
