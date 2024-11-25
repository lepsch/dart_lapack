// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zla_wwaddw(
  final int N,
  final Array<Complex> X_,
  final Array<Complex> Y_,
  final Array<Complex> W_,
) {
  final X = X_.having();
  final Y = Y_.having();
  final W = W_.having();

  for (var I = 1; I <= N; I++) {
    var S = X[I] + W[I];
    S = (S + S) - S;
    Y[I] = ((X[I] - S) + W[I]) + Y[I];
    X[I] = S;
  }
}
