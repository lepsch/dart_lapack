// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dsyr2k(
  final String UPLO,
  final String TRANS,
  final int N,
  final int K,
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

  // Test the input parameters.

  final NROWA = lsame(TRANS, 'N') ? N : K;
  final UPPER = lsame(UPLO, 'U');

  var INFO = 0;
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO = 2;
  } else if (N < 0) {
    INFO = 3;
  } else if (K < 0) {
    INFO = 4;
  } else if (LDA < max(1, NROWA)) {
    INFO = 7;
  } else if (LDB < max(1, NROWA)) {
    INFO = 9;
  } else if (LDC < max(1, N)) {
    INFO = 12;
  }
  if (INFO != 0) {
    xerbla('DSYR2K', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) return;

  // And when  alpha == zero.

  if (ALPHA == ZERO) {
    if (UPPER) {
      if (BETA == ZERO) {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= J; I++) {
            C[I][J] = ZERO;
          }
        }
      } else {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= J; I++) {
            C[I][J] *= BETA;
          }
        }
      }
    } else {
      if (BETA == ZERO) {
        for (var J = 1; J <= N; J++) {
          for (var I = J; I <= N; I++) {
            C[I][J] = ZERO;
          }
        }
      } else {
        for (var J = 1; J <= N; J++) {
          for (var I = J; I <= N; I++) {
            C[I][J] *= BETA;
          }
        }
      }
    }
    return;
  }

  // Start the operations.

  if (lsame(TRANS, 'N')) {
    // Form  C := alpha*A*B**T + alpha*B*A**T + C.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (var I = 1; I <= J; I++) {
            C[I][J] = ZERO;
          }
        } else if (BETA != ONE) {
          for (var I = 1; I <= J; I++) {
            C[I][J] *= BETA;
          }
        }
        for (var L = 1; L <= K; L++) {
          if ((A[J][L] != ZERO) || (B[J][L] != ZERO)) {
            final TEMP1 = ALPHA * B[J][L];
            final TEMP2 = ALPHA * A[J][L];
            for (var I = 1; I <= J; I++) {
              C[I][J] += A[I][L] * TEMP1 + B[I][L] * TEMP2;
            }
          }
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (var I = J; I <= N; I++) {
            C[I][J] = ZERO;
          }
        } else if (BETA != ONE) {
          for (var I = J; I <= N; I++) {
            C[I][J] *= BETA;
          }
        }
        for (var L = 1; L <= K; L++) {
          if ((A[J][L] != ZERO) || (B[J][L] != ZERO)) {
            final TEMP1 = ALPHA * B[J][L];
            final TEMP2 = ALPHA * A[J][L];
            for (var I = J; I <= N; I++) {
              C[I][J] += A[I][L] * TEMP1 + B[I][L] * TEMP2;
            }
          }
        }
      }
    }
  } else {
    // Form  C := alpha*A**T*B + alpha*B**T*A + C.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J; I++) {
          var TEMP1 = ZERO, TEMP2 = ZERO;
          for (var L = 1; L <= K; L++) {
            TEMP1 += A[L][I] * B[L][J];
            TEMP2 += B[L][I] * A[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA * TEMP1 + ALPHA * TEMP2;
          } else {
            C[I][J] = BETA * C[I][J] + ALPHA * TEMP1 + ALPHA * TEMP2;
          }
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = J; I <= N; I++) {
          var TEMP1 = ZERO, TEMP2 = ZERO;
          for (var L = 1; L <= K; L++) {
            TEMP1 += A[L][I] * B[L][J];
            TEMP2 += B[L][I] * A[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA * TEMP1 + ALPHA * TEMP2;
          } else {
            C[I][J] = BETA * C[I][J] + ALPHA * TEMP1 + ALPHA * TEMP2;
          }
        }
      }
    }
  }
}
