// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/matrix.dart';

void dgtts2(
  final int ITRANS,
  final int N,
  final int NRHS,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final B = B_.having(ld: LDB);
  int I, IP, J;
  double TEMP;

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (ITRANS == 0) {
    // Solve A*X = B using the LU factorization of A,
    // overwriting each right hand side vector with its solution.

    if (NRHS <= 1) {
      for (J = 1; J <= NRHS; J++) {
        // Solve L*x = b.

        for (I = 1; I <= N - 1; I++) {
          IP = IPIV[I];
          TEMP = B[I + 1 - IP + I][J] - DL[I] * B[IP][J];
          B[I][J] = B[IP][J];
          B[I + 1][J] = TEMP;
        }

        // Solve U*x = b.

        B[N][J] /= D[N];
        if (N > 1) B[N - 1][J] = (B[N - 1][J] - DU[N - 1] * B[N][J]) / D[N - 1];
        for (I = N - 2; I >= 1; I--) {
          B[I][J] =
              (B[I][J] - DU[I] * B[I + 1][J] - DU2[I] * B[I + 2][J]) / D[I];
        }
      }
    } else {
      for (J = 1; J <= NRHS; J++) {
        // Solve L*x = b.

        for (I = 1; I <= N - 1; I++) {
          if (IPIV[I] == I) {
            B[I + 1][J] -= DL[I] * B[I][J];
          } else {
            TEMP = B[I][J];
            B[I][J] = B[I + 1][J];
            B[I + 1][J] = TEMP - DL[I] * B[I][J];
          }
        }

        // Solve U*x = b.

        B[N][J] /= D[N];
        if (N > 1) B[N - 1][J] = (B[N - 1][J] - DU[N - 1] * B[N][J]) / D[N - 1];
        for (I = N - 2; I >= 1; I--) {
          B[I][J] =
              (B[I][J] - DU[I] * B[I + 1][J] - DU2[I] * B[I + 2][J]) / D[I];
        }
      }
    }
  } else {
    // Solve A**T * X = B.

    if (NRHS <= 1) {
      // Solve U**T*x = b.

      for (J = 1; J <= NRHS; J++) {
        B[1][J] /= D[1];
        if (N > 1) B[2][J] = (B[2][J] - DU[1] * B[1][J]) / D[2];
        for (I = 3; I <= N; I++) {
          B[I][J] =
              (B[I][J] - DU[I - 1] * B[I - 1][J] - DU2[I - 2] * B[I - 2][J]) /
                  D[I];
        }

        // Solve L**T*x = b.

        for (I = N - 1; I >= 1; I--) {
          IP = IPIV[I];
          TEMP = B[I][J] - DL[I] * B[I + 1][J];
          B[I][J] = B[IP][J];
          B[IP][J] = TEMP;
        }
      }
    } else {
      for (J = 1; J <= NRHS; J++) {
        // Solve U**T*x = b.

        B[1][J] /= D[1];
        if (N > 1) B[2][J] = (B[2][J] - DU[1] * B[1][J]) / D[2];
        for (I = 3; I <= N; I++) {
          B[I][J] =
              (B[I][J] - DU[I - 1] * B[I - 1][J] - DU2[I - 2] * B[I - 2][J]) /
                  D[I];
        }
        for (I = N - 1; I >= 1; I--) {
          if (IPIV[I] == I) {
            B[I][J] -= DL[I] * B[I + 1][J];
          } else {
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - DL[I] * TEMP;
            B[I][J] = TEMP;
          }
        }
      }
    }
  }
}
