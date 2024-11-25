// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlaptm(
  final int N,
  final int NRHS,
  final double ALPHA,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> X_,
  final int LDX,
  final double BETA,
  final Matrix<double> B_,
  final int LDB,
) {
  final D = D_.having();
  final E = E_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  const ONE = 1.0, ZERO = 0.0;

  if (N == 0) return;

  // Multiply B by BETA if BETA != 1.

  if (BETA == ZERO) {
    for (var J = 1; J <= NRHS; J++) {
      for (var I = 1; I <= N; I++) {
        B[I][J] = ZERO;
      }
    }
  } else if (BETA == -ONE) {
    for (var J = 1; J <= NRHS; J++) {
      for (var I = 1; I <= N; I++) {
        B[I][J] = -B[I][J];
      }
    }
  }

  if (ALPHA == ONE) {
    // Compute B := B + A*X

    for (var J = 1; J <= NRHS; J++) {
      if (N == 1) {
        B[1][J] += D[1] * X[1][J];
      } else {
        B[1][J] += D[1] * X[1][J] + E[1] * X[2][J];
        B[N][J] += E[N - 1] * X[N - 1][J] + D[N] * X[N][J];
        for (var I = 2; I <= N - 1; I++) {
          B[I][J] +=
              E[I - 1] * X[I - 1][J] + D[I] * X[I][J] + E[I] * X[I + 1][J];
        }
      }
    }
  } else if (ALPHA == -ONE) {
    // Compute B := B - A*X

    for (var J = 1; J <= NRHS; J++) {
      if (N == 1) {
        B[1][J] -= D[1] * X[1][J];
      } else {
        B[1][J] -= D[1] * X[1][J] + E[1] * X[2][J];
        B[N][J] -= E[N - 1] * X[N - 1][J] + D[N] * X[N][J];
        for (var I = 2; I <= N - 1; I++) {
          B[I][J] -=
              E[I - 1] * X[I - 1][J] + D[I] * X[I][J] + E[I] * X[I + 1][J];
        }
      }
    }
  }
}
