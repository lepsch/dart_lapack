// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/matrix.dart';

void dptts2(
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> B_,
  final int LDB,
) {
  final D = D_.having();
  final E = E_.having();
  final B = B_.having(ld: LDB);
  int I, J;

  // Quick return if possible

  if (N <= 1) {
    if (N == 1) dscal(NRHS, 1.0 / D[1], B.asArray(), LDB);
    return;
  }

  // Solve A * X = B using the factorization A = L*D*L**T,
  // overwriting each right hand side vector with its solution.

  for (J = 1; J <= NRHS; J++) {
    // Solve L * x = b.

    for (I = 2; I <= N; I++) {
      B[I][J] -= B[I - 1][J] * E[I - 1];
    }

    // Solve D * L**T * x = b.

    B[N][J] /= D[N];
    for (I = N - 1; I >= 1; I--) {
      B[I][J] = B[I][J] / D[I] - B[I + 1][J] * E[I];
    }
  }
}
