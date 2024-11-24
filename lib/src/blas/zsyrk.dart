// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zsyrk(
  final String UPLO,
  final String TRANS,
  final int N,
  final int K,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Complex BETA,
  final Matrix<Complex> C_,
  final int LDC,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final C = C_.having(ld: LDC);
  Complex TEMP;
  int I, INFO, J, L, NROWA;
  bool UPPER;

  // Test the input parameters.

  if (lsame(TRANS, 'N')) {
    NROWA = N;
  } else {
    NROWA = K;
  }
  UPPER = lsame(UPLO, 'U');

  INFO = 0;
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T')) {
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
    xerbla('ZSYRK', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) ||
      (((ALPHA == Complex.zero) || (K == 0)) && (BETA == Complex.one))) return;

  // And when  alpha == zero.

  if (ALPHA == Complex.zero) {
    if (UPPER) {
      if (BETA == Complex.zero) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            C[I][J] = Complex.zero;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
      }
    } else {
      if (BETA == Complex.zero) {
        for (J = 1; J <= N; J++) {
          for (I = J; I <= N; I++) {
            C[I][J] = Complex.zero;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = J; I <= N; I++) {
            C[I][J] = BETA * C[I][J];
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
      for (J = 1; J <= N; J++) {
        if (BETA == Complex.zero) {
          for (I = 1; I <= J; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != Complex.one) {
          for (I = 1; I <= J; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          if (A[J][L] != Complex.zero) {
            TEMP = ALPHA * A[J][L];
            for (I = 1; I <= J; I++) {
              C[I][J] += TEMP * A[I][L];
            }
          }
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        if (BETA == Complex.zero) {
          for (I = J; I <= N; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != Complex.one) {
          for (I = J; I <= N; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          if (A[J][L] != Complex.zero) {
            TEMP = ALPHA * A[J][L];
            for (I = J; I <= N; I++) {
              C[I][J] += TEMP * A[I][L];
            }
          }
        }
      }
    }
  } else {
    // Form  C := alpha*A**T*A + beta*C.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I] * A[L][J];
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        for (I = J; I <= N; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I] * A[L][J];
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    }
  }
}
