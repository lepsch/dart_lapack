// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlartv(
  final int N,
  final Array<Complex> X_,
  final int INCX,
  final Array<Complex> Y_,
  final int INCY,
  final Array<double> C_,
  final Array<Complex> S_,
  final int INCC,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final Y = Y_.having();
  final C = C_.having();
  final S = S_.having();
  int I, IC, IX, IY;
  Complex XI, YI;

  IX = 1;
  IY = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    XI = X[IX];
    YI = Y[IY];
    X[IX] = C[IC].toComplex() * XI + S[IC] * YI;
    Y[IY] = C[IC].toComplex() * YI - S[IC].conjugate() * XI;
    IX += INCX;
    IY += INCY;
    IC += INCC;
  }
}
