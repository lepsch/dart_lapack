// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlaruv.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlarnv(
  final int IDIST,
  final Array<int> ISEED_,
  final int N,
  final Array<Complex> X_,
) {
  final ISEED = ISEED_.having(length: 4);
  final X = X_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const LV = 128;
  const TWOPI = 6.28318530717958647692528676655900576839;
  int I, IL, IV;
  final U = Array<double>(LV);

  for (IV = 1; IV <= N; IV += LV ~/ 2) {
    IL = min(LV ~/ 2, N - IV + 1);

    // Call DLARUV to generate 2*IL real numbers from a uniform (0,1)
    // distribution (2*IL <= LV)

    dlaruv(ISEED, 2 * IL, U);

    if (IDIST == 1) {
      // Copy generated numbers

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = Complex(U[2 * I - 1], U[2 * I]);
      }
    } else if (IDIST == 2) {
      // Convert generated numbers to uniform (-1,1) distribution

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = Complex(TWO * U[2 * I - 1] - ONE, TWO * U[2 * I] - ONE);
      }
    } else if (IDIST == 3) {
      // Convert generated numbers to normal (0,1) distribution

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = sqrt(-TWO * log(U[2 * I - 1])).toComplex() *
            Complex(ZERO, TWOPI * U[2 * I]).exp();
      }
    } else if (IDIST == 4) {
      // Convert generated numbers to complex numbers uniformly
      // distributed on the unit disk

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = sqrt(U[2 * I - 1]).toComplex() *
            Complex(ZERO, TWOPI * U[2 * I]).exp();
      }
    } else if (IDIST == 5) {
      // Convert generated numbers to complex numbers uniformly
      // distributed on the unit circle

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = Complex(ZERO, TWOPI * U[2 * I]).exp();
      }
    }
  }
}
