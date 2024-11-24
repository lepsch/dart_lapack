// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dsymm(
  final String SIDE,
  final String UPLO,
  final int M,
  final int N,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final double BETA,
  final Matrix<double> C_,
  final int LDC,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  const ONE = 1.0, ZERO = 0.0;

  // Set NROWA as the number of rows of A.

  final NROWA = lsame(SIDE, 'L') ? M : N;
  final UPPER = lsame(UPLO, 'U');

  // Test the input parameters.

  var INFO = 0;
  if (!lsame(SIDE, 'L') && !lsame(SIDE, 'R')) {
    INFO = 1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO = 2;
  } else if (M < 0) {
    INFO = 3;
  } else if (N < 0) {
    INFO = 4;
  } else if (LDA < max(1, NROWA)) {
    INFO = 7;
  } else if (LDB < max(1, M)) {
    INFO = 9;
  } else if (LDC < max(1, M)) {
    INFO = 12;
  }
  if (INFO != 0) {
    xerbla('DSYMM', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // And when  alpha == zero.

  if (ALPHA == ZERO) {
    if (BETA == ZERO) {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= M; I++) {
          C[I][J] = ZERO;
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= M; I++) {
          C[I][J] *= BETA;
        }
      }
    }
    return;
  }

  // Start the operations.

  if (lsame(SIDE, 'L')) {
    // Form  C := alpha*A*B + beta*C.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= M; I++) {
          final TEMP1 = ALPHA * B[I][J];
          var TEMP2 = ZERO;
          for (var K = 1; K <= I - 1; K++) {
            C[K][J] += TEMP1 * A[K][I];
            TEMP2 += B[K][J] * A[K][I];
          }
          if (BETA == ZERO) {
            C[I][J] = TEMP1 * A[I][I] + ALPHA * TEMP2;
          } else {
            C[I][J] = BETA * C[I][J] + TEMP1 * A[I][I] + ALPHA * TEMP2;
          }
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = M; I >= 1; I--) {
          final TEMP1 = ALPHA * B[I][J];
          var TEMP2 = ZERO;
          for (var K = I + 1; K <= M; K++) {
            C[K][J] += TEMP1 * A[K][I];
            TEMP2 += B[K][J] * A[K][I];
          }
          if (BETA == ZERO) {
            C[I][J] = TEMP1 * A[I][I] + ALPHA * TEMP2;
          } else {
            C[I][J] = BETA * C[I][J] + TEMP1 * A[I][I] + ALPHA * TEMP2;
          }
        }
      }
    }
  } else {
    // Form  C := alpha*B*A + beta*C.

    for (var J = 1; J <= N; J++) {
      final TEMP1 = ALPHA * A[J][J];
      if (BETA == ZERO) {
        for (var I = 1; I <= M; I++) {
          C[I][J] = TEMP1 * B[I][J];
        }
      } else {
        for (var I = 1; I <= M; I++) {
          C[I][J] = BETA * C[I][J] + TEMP1 * B[I][J];
        }
      }
      for (var K = 1; K <= J - 1; K++) {
        final double TEMP1;
        if (UPPER) {
          TEMP1 = ALPHA * A[K][J];
        } else {
          TEMP1 = ALPHA * A[J][K];
        }
        for (var I = 1; I <= M; I++) {
          C[I][J] += TEMP1 * B[I][K];
        }
      }
      for (var K = J + 1; K <= N; K++) {
        final double TEMP1;
        if (UPPER) {
          TEMP1 = ALPHA * A[J][K];
        } else {
          TEMP1 = ALPHA * A[K][J];
        }
        for (var I = 1; I <= M; I++) {
          C[I][J] += TEMP1 * B[I][K];
        }
      }
    }
  }
}
