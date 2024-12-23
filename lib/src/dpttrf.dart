// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dpttrf(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0;
  int I, I4;
  double EI;

  // Test the input parameters.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('DPTTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Compute the L*D*L**T (or U**T*D*U) factorization of A.

  I4 = ((N - 1) % 4);
  for (I = 1; I <= I4; I++) {
    if (D[I] <= ZERO) {
      INFO.value = I;
      return;
    }
    EI = E[I];
    E[I] = EI / D[I];
    D[I + 1] -= E[I] * EI;
  }

  for (I = I4 + 1; I <= N - 4; I += 4) {
    // Drop out of the loop if d[i] <= 0: the matrix is not positive
    // definite.

    if (D[I] <= ZERO) {
      INFO.value = I;
      return;
    }

    // Solve for e[i] and d[i+1].

    EI = E[I];
    E[I] = EI / D[I];
    D[I + 1] -= E[I] * EI;

    if (D[I + 1] <= ZERO) {
      INFO.value = I + 1;
      return;
    }

    // Solve for e[i+1] and d[i+2].

    EI = E[I + 1];
    E[I + 1] = EI / D[I + 1];
    D[I + 2] -= E[I + 1] * EI;

    if (D[I + 2] <= ZERO) {
      INFO.value = I + 2;
      return;
    }

    // Solve for e[i+2] and d[i+3].

    EI = E[I + 2];
    E[I + 2] = EI / D[I + 2];
    D[I + 3] -= E[I + 2] * EI;

    if (D[I + 3] <= ZERO) {
      INFO.value = I + 3;
      return;
    }

    // Solve for e[i+3] and d[i+4].

    EI = E[I + 3];
    E[I + 3] = EI / D[I + 3];
    D[I + 4] -= E[I + 3] * EI;
  }

  // Check d[n] for positive definiteness.

  if (D[N] <= ZERO) INFO.value = N;
}
