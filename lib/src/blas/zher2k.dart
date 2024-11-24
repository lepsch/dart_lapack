// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zher2k(
  final String UPLO,
  final String TRANS,
  final int N,
  final int K,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final double BETA,
  final Matrix<Complex> C_,
  final int LDC,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  Complex TEMP1, TEMP2;
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
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'C')) {
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
    xerbla('ZHER2K', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (((ALPHA == Complex.zero) || (K == 0)) && (BETA == ONE))) {
    return;
  }

  // And when  alpha == zero.

  if (ALPHA == Complex.zero) {
    if (UPPER) {
      if (BETA == ZERO) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            C[I][J] = Complex.zero;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J - 1; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
          C[J][J] = (BETA * C[J][J].real).toComplex();
        }
      }
    } else {
      if (BETA == ZERO) {
        for (J = 1; J <= N; J++) {
          for (I = J; I <= N; I++) {
            C[I][J] = Complex.zero;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          C[J][J] = (BETA * C[J][J].real).toComplex();
          for (I = J + 1; I <= N; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
        }
      }
    }
    return;
  }

  // Start the operations.

  if (lsame(TRANS, 'N')) {
    // Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
    // C.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (I = 1; I <= J; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != ONE) {
          for (I = 1; I <= J - 1; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
          C[J][J] = (BETA * C[J][J].real).toComplex();
        } else {
          C[J][J] = C[J][J].real.toComplex();
        }
        for (L = 1; L <= K; L++) {
          if ((A[J][L] != Complex.zero) || (B[J][L] != Complex.zero)) {
            TEMP1 = ALPHA * B[J][L].conjugate();
            TEMP2 = (ALPHA * A[J][L]).conjugate();
            for (I = 1; I <= J - 1; I++) {
              C[I][J] += A[I][L] * TEMP1 + B[I][L] * TEMP2;
            }
            C[J][J] = (C[J][J].real + (A[J][L] * TEMP1 + B[J][L] * TEMP2).real)
                .toComplex();
          }
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (I = J; I <= N; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != ONE) {
          for (I = J + 1; I <= N; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
          C[J][J] = (BETA * C[J][J].real).toComplex();
        } else {
          C[J][J] = C[J][J].real.toComplex();
        }
        for (L = 1; L <= K; L++) {
          if ((A[J][L] != Complex.zero) || (B[J][L] != Complex.zero)) {
            TEMP1 = ALPHA * B[J][L].conjugate();
            TEMP2 = (ALPHA * A[J][L]).conjugate();
            for (I = J + 1; I <= N; I++) {
              C[I][J] += A[I][L] * TEMP1 + B[I][L] * TEMP2;
            }
            C[J][J] = (C[J][J].real + (A[J][L] * TEMP1 + B[J][L] * TEMP2).real)
                .toComplex();
          }
        }
      }
    }
  } else {
    // Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
    // C.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          TEMP1 = Complex.zero;
          TEMP2 = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP1 += A[L][I].conjugate() * B[L][J];
            TEMP2 += B[L][I].conjugate() * A[L][J];
          }
          if (I == J) {
            if (BETA == ZERO) {
              C[J][J] =
                  (ALPHA * TEMP1 + ALPHA.conjugate() * TEMP2).real.toComplex();
            } else {
              C[J][J] = (BETA * C[J][J].real +
                      (ALPHA * TEMP1 + ALPHA.conjugate() * TEMP2).real)
                  .toComplex();
            }
          } else {
            if (BETA == ZERO) {
              C[I][J] = ALPHA * TEMP1 + ALPHA.conjugate() * TEMP2;
            } else {
              C[I][J] = BETA.toComplex() * C[I][J] +
                  ALPHA * TEMP1 +
                  ALPHA.conjugate() * TEMP2;
            }
          }
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        for (I = J; I <= N; I++) {
          TEMP1 = Complex.zero;
          TEMP2 = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP1 += A[L][I].conjugate() * B[L][J];
            TEMP2 += B[L][I].conjugate() * A[L][J];
          }
          if (I == J) {
            if (BETA == ZERO) {
              C[J][J] =
                  (ALPHA * TEMP1 + ALPHA.conjugate() * TEMP2).real.toComplex();
            } else {
              C[J][J] = (BETA * C[J][J].real +
                      (ALPHA * TEMP1 + ALPHA.conjugate() * TEMP2).real)
                  .toComplex();
            }
          } else {
            if (BETA == ZERO) {
              C[I][J] = ALPHA * TEMP1 + ALPHA.conjugate() * TEMP2;
            } else {
              C[I][J] = BETA.toComplex() * C[I][J] +
                  ALPHA * TEMP1 +
                  ALPHA.conjugate() * TEMP2;
            }
          }
        }
      }
    }
  }
}
