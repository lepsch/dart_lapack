// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlascl2(
  final int M,
  final int N,
  final Array<double> D,
  final Matrix<Complex> X_,
  final int LDX,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having(ld: LDX);

  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= M; I++) {
      X[I][J] *= D[I].toComplex();
    }
  }
}
