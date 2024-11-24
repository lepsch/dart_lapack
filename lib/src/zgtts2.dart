// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zgtts2(
  final int ITRANS,
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DU2_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  int I, J = 0;
  Complex TEMP;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC DCONJG

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (ITRANS == 0) {
    // Solve A*X = B using the LU factorization of A,
    // overwriting each right hand side vector with its solution.

    if (NRHS <= 1) {
      J = 1;
      while (true) {
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
        if (J < NRHS) {
          J++;
          continue;
        }
        break;
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
  } else if (ITRANS == 1) {
    // Solve A**T * X = B.

    if (NRHS <= 1) {
      J = 1;
      while (true) {
        // Solve U**T * x = b.

        B[1][J] /= D[1];
        if (N > 1) B[2][J] = (B[2][J] - DU[1] * B[1][J]) / D[2];
        for (I = 3; I <= N; I++) {
          B[I][J] =
              (B[I][J] - DU[I - 1] * B[I - 1][J] - DU2[I - 2] * B[I - 2][J]) /
                  D[I];
        }

        // Solve L**T * x = b.

        for (I = N - 1; I >= 1; I--) {
          if (IPIV[I] == I) {
            B[I][J] -= DL[I] * B[I + 1][J];
          } else {
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - DL[I] * TEMP;
            B[I][J] = TEMP;
          }
        }
        if (J < NRHS) {
          J++;
          continue;
        }
        break;
      }
    } else {
      for (J = 1; J <= NRHS; J++) {
        // Solve U**T * x = b.

        B[1][J] /= D[1];
        if (N > 1) B[2][J] = (B[2][J] - DU[1] * B[1][J]) / D[2];
        for (I = 3; I <= N; I++) {
          B[I][J] =
              (B[I][J] - DU[I - 1] * B[I - 1][J] - DU2[I - 2] * B[I - 2][J]) /
                  D[I];
        }

        // Solve L**T * x = b.

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
  } else {
    // Solve A**H * X = B.

    if (NRHS <= 1) {
      J = 1;
      while (true) {
        // Solve U**H * x = b.

        B[1][J] /= D[1].conjugate();
        if (N > 1) {
          B[2][J] = (B[2][J] - DU[1].conjugate() * B[1][J]) / D[2].conjugate();
        }
        for (I = 3; I <= N; I++) {
          B[I][J] = (B[I][J] -
                  DU[I - 1].conjugate() * B[I - 1][J] -
                  DU2[I - 2].conjugate() * B[I - 2][J]) /
              D[I].conjugate();
        }

        // Solve L**H * x = b.

        for (I = N - 1; I >= 1; I--) {
          if (IPIV[I] == I) {
            B[I][J] -= DL[I].conjugate() * B[I + 1][J];
          } else {
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - DL[I].conjugate() * TEMP;
            B[I][J] = TEMP;
          }
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

        B[1][J] /= D[1].conjugate();
        if (N > 1) {
          B[2][J] = (B[2][J] - DU[1].conjugate() * B[1][J]) / D[2].conjugate();
        }
        for (I = 3; I <= N; I++) {
          B[I][J] = (B[I][J] -
                  DU[I - 1].conjugate() * B[I - 1][J] -
                  DU2[I - 2].conjugate() * B[I - 2][J]) /
              D[I].conjugate();
        }

        // Solve L**H * x = b.

        for (I = N - 1; I >= 1; I--) {
          if (IPIV[I] == I) {
            B[I][J] -= DL[I].conjugate() * B[I + 1][J];
          } else {
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - DL[I].conjugate() * TEMP;
            B[I][J] = TEMP;
          }
        }
      }
    }
  }
}
