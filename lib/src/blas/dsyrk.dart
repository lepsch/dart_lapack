// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dsyrk(
  final String UPLO,
  final String TRANS,
  final int N,
  final int K,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final double BETA,
  final Matrix<double> C_,
  final int LDC,
) {
  final A = A_.having(ld: LDA);
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
  } else if (LDC < max(1, N)) {
    INFO = 10;
  }
  if (INFO != 0) {
    xerbla('DSYRK', INFO);
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
    // Form  C := alpha*A*A**T + beta*C.

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
          if (A[J][L] != ZERO) {
            final TEMP = ALPHA * A[J][L];
            for (var I = 1; I <= J; I++) {
              C[I][J] += TEMP * A[I][L];
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
          if (A[J][L] != ZERO) {
            final TEMP = ALPHA * A[J][L];
            for (var I = J; I <= N; I++) {
              C[I][J] += TEMP * A[I][L];
            }
          }
        }
      }
    }
  } else {
    // Form  C := alpha*A**T*A + beta*C.

    if (UPPER) {
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= J; I++) {
          var TEMP = ZERO;
          for (var L = 1; L <= K; L++) {
            TEMP += A[L][I] * A[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        for (var I = J; I <= N; I++) {
          var TEMP = ZERO;
          for (var L = 1; L <= K; L++) {
            TEMP += A[L][I] * A[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    }
  }
}
