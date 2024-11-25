// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/dlaruv.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarnv(
  final int IDIST,
  final Array<int> ISEED_,
  final int N,
  final Array<double> X_,
) {
  final X = X_.having();
  final ISEED = ISEED_.having();
  const ONE = 1.0, TWO = 2.0;
  const LV = 128;
  const TWOPI = 6.28318530717958647692528676655900576839;
  int I, IL, IL2, IV;
  final U = Array<double>(LV);

  for (IV = 1; IV <= N; IV += LV ~/ 2) {
    IL = min(LV ~/ 2, N - IV + 1);
    if (IDIST == 3) {
      IL2 = 2 * IL;
    } else {
      IL2 = IL;
    }

    // Call DLARUV to generate IL2 numbers from a uniform (0,1)
    // distribution (IL2 <= LV)

    dlaruv(ISEED, IL2, U);

    if (IDIST == 1) {
      // Copy generated numbers

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = U[I];
      }
    } else if (IDIST == 2) {
      // Convert generated numbers to uniform (-1,1) distribution

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = TWO * U[I] - ONE;
      }
    } else if (IDIST == 3) {
      // Convert generated numbers to normal (0,1) distribution

      for (I = 1; I <= IL; I++) {
        X[IV + I - 1] = sqrt(-TWO * log(U[2 * I - 1])) * cos(TWOPI * U[2 * I]);
      }
    }
  }
}
