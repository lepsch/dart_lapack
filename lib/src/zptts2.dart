// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zptts2(
  final int IUPLO,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.having(ld: LDB);
  final E = E_.having();
  final D = D_.having();
  int I, J;

  // Quick return if possible

  if (N <= 1) {
    if (N == 1) zdscal(NRHS, 1.0 / D[1], B.asArray(), LDB);
    return;
  }

  if (IUPLO == 1) {
    // Solve A * X = B using the factorization A = U**H *D*U,
    // overwriting each right hand side vector with its solution.

    if (NRHS <= 2) {
      J = 1;
      while (true) {
        // Solve U**H * x = b.

        for (I = 2; I <= N; I++) {
          B[I][J] -= B[I - 1][J] * E[I - 1].conjugate();
        }

        // Solve D * U * x = b.

        for (I = 1; I <= N; I++) {
          B[I][J] /= D[I].toComplex();
        }
        for (I = N - 1; I >= 1; I--) {
          B[I][J] -= B[I + 1][J] * E[I];
        }
        if (J < NRHS) {
          J++;
          continue;
        }
        break;
      }
    } else {
      for (J = 1; J <= NRHS; J++) {
        // Solve U**H * x = b.

        for (I = 2; I <= N; I++) {
          B[I][J] -= B[I - 1][J] * E[I - 1].conjugate();
        }

        // Solve D * U * x = b.

        B[N][J] /= D[N].toComplex();
        for (I = N - 1; I >= 1; I--) {
          B[I][J] = B[I][J] / D[I].toComplex() - B[I + 1][J] * E[I];
        }
      }
    }
  } else {
    // Solve A * X = B using the factorization A = L*D*L**H,
    // overwriting each right hand side vector with its solution.

    if (NRHS <= 2) {
      J = 1;
      while (true) {
        // Solve L * x = b.

        for (I = 2; I <= N; I++) {
          B[I][J] -= B[I - 1][J] * E[I - 1];
        }

        // Solve D * L**H * x = b.

        for (I = 1; I <= N; I++) {
          B[I][J] /= D[I].toComplex();
        }
        for (I = N - 1; I >= 1; I--) {
          B[I][J] -= B[I + 1][J] * E[I].conjugate();
        }
        if (J < NRHS) {
          J++;
          continue;
        }
        break;
      }
    } else {
      for (J = 1; J <= NRHS; J++) {
        // Solve L * x = b.

        for (I = 2; I <= N; I++) {
          B[I][J] -= B[I - 1][J] * E[I - 1];
        }

        // Solve D * L**H * x = b.

        B[N][J] /= D[N].toComplex();
        for (I = N - 1; I >= 1; I--) {
          B[I][J] = B[I][J] / D[I].toComplex() - B[I + 1][J] * E[I].conjugate();
        }
      }
    }
  }
}
