// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zpttrf(
  final int N,
  final Array<double> D_,
  final Array<Complex> E_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0;
  int I, I4;
  double EII, EIR, F, G;

  // Test the input parameters.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('ZPTTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Compute the L*D*L**H (or U**H *D*U) factorization of A.

  I4 = ((N - 1) % 4);
  for (I = 1; I <= I4; I++) {
    if (D[I] <= ZERO) {
      INFO.value = I;
      return;
    }
    EIR = E[I].real;
    EII = E[I].imaginary;
    F = EIR / D[I];
    G = EII / D[I];
    E[I] = Complex(F, G);
    D[I + 1] -= F * EIR + G * EII;
  }

  for (I = I4 + 1; I <= N - 4; I += 4) {
    // Drop out of the loop if d(i) <= 0: the matrix is not positive
    // definite.

    if (D[I] <= ZERO) {
      INFO.value = I;
      return;
    }

    // Solve for e(i) and d(i+1).

    EIR = E[I].real;
    EII = E[I].imaginary;
    F = EIR / D[I];
    G = EII / D[I];
    E[I] = Complex(F, G);
    D[I + 1] -= F * EIR + G * EII;

    if (D[I + 1] <= ZERO) {
      INFO.value = I + 1;
      return;
    }

    // Solve for e(i+1) and d(i+2).

    EIR = E[I + 1].real;
    EII = E[I + 1].imaginary;
    F = EIR / D[I + 1];
    G = EII / D[I + 1];
    E[I + 1] = Complex(F, G);
    D[I + 2] -= F * EIR + G * EII;

    if (D[I + 2] <= ZERO) {
      INFO.value = I + 2;
      return;
    }

    // Solve for e(i+2) and d(i+3).

    EIR = E[I + 2].real;
    EII = E[I + 2].imaginary;
    F = EIR / D[I + 2];
    G = EII / D[I + 2];
    E[I + 2] = Complex(F, G);
    D[I + 3] -= F * EIR + G * EII;

    if (D[I + 3] <= ZERO) {
      INFO.value = I + 3;
      return;
    }

    // Solve for e(i+3) and d(i+4).

    EIR = E[I + 3].real;
    EII = E[I + 3].imaginary;
    F = EIR / D[I + 3];
    G = EII / D[I + 3];
    E[I + 3] = Complex(F, G);
    D[I + 4] -= F * EIR + G * EII;
  }

  // Check d(n) for positive definiteness.

  if (D[N] <= ZERO) INFO.value = N;
}
