// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zgemm(
  final String TRANSA,
  final String TRANSB,
  final int M,
  final int N,
  final int K,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Complex BETA,
  final Matrix<Complex> C_,
  final int LDC,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  Complex TEMP;
  int I, INFO, J, L, NROWA, NROWB;
  bool CONJA, CONJB, NOTA, NOTB;

  // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  // conjugated or transposed, set  CONJA and CONJB  as true if  A  and
  // B  respectively are to be  transposed but  not conjugated  and set
  // NROWA and NROWB  as the number of rows  of  A  and  B  respectively.

  NOTA = lsame(TRANSA, 'N');
  NOTB = lsame(TRANSB, 'N');
  CONJA = lsame(TRANSA, 'C');
  CONJB = lsame(TRANSB, 'C');
  if (NOTA) {
    NROWA = M;
  } else {
    NROWA = K;
  }
  if (NOTB) {
    NROWB = K;
  } else {
    NROWB = N;
  }

  // Test the input parameters.

  INFO = 0;
  if (!NOTA && !CONJA && !lsame(TRANSA, 'T')) {
    INFO = 1;
  } else if (!NOTB && !CONJB && !lsame(TRANSB, 'T')) {
    INFO = 2;
  } else if (M < 0) {
    INFO = 3;
  } else if (N < 0) {
    INFO = 4;
  } else if (K < 0) {
    INFO = 5;
  } else if (LDA < max(1, NROWA)) {
    INFO = 8;
  } else if (LDB < max(1, NROWB)) {
    INFO = 10;
  } else if (LDC < max(1, M)) {
    INFO = 13;
  }
  if (INFO != 0) {
    xerbla('ZGEMM', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) ||
      (N == 0) ||
      (((ALPHA == Complex.zero) || (K == 0)) && (BETA == Complex.one))) return;

  // And when  alpha == zero.

  if (ALPHA == Complex.zero) {
    if (BETA == Complex.zero) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          C[I][J] = Complex.zero;
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          C[I][J] = BETA * C[I][J];
        }
      }
    }
    return;
  }

  // Start the operations.

  if (NOTB) {
    if (NOTA) {
      // Form  C := alpha*A*B + beta*C.

      for (J = 1; J <= N; J++) {
        if (BETA == Complex.zero) {
          for (I = 1; I <= M; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != Complex.one) {
          for (I = 1; I <= M; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          TEMP = ALPHA * B[L][J];
          for (I = 1; I <= M; I++) {
            C[I][J] += TEMP * A[I][L];
          }
        }
      }
    } else if (CONJA) {
      // Form  C := alpha*A**H*B + beta*C.

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I].conjugate() * B[L][J];
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    } else {
      // Form  C := alpha*A**T*B + beta*C

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I] * B[L][J];
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    }
  } else if (NOTA) {
    if (CONJB) {
      // Form  C := alpha*A*B**H + beta*C.

      for (J = 1; J <= N; J++) {
        if (BETA == Complex.zero) {
          for (I = 1; I <= M; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != Complex.one) {
          for (I = 1; I <= M; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          TEMP = ALPHA * B[J][L].conjugate();
          for (I = 1; I <= M; I++) {
            C[I][J] += TEMP * A[I][L];
          }
        }
      }
    } else {
      // Form  C := alpha*A*B**T + beta*C

      for (J = 1; J <= N; J++) {
        if (BETA == Complex.zero) {
          for (I = 1; I <= M; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != Complex.one) {
          for (I = 1; I <= M; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          TEMP = ALPHA * B[J][L];
          for (I = 1; I <= M; I++) {
            C[I][J] += TEMP * A[I][L];
          }
        }
      }
    }
  } else if (CONJA) {
    if (CONJB) {
      // Form  C := alpha*A**H*B**H + beta*C.

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I].conjugate() * B[J][L].conjugate();
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    } else {
      // Form  C := alpha*A**H*B**T + beta*C

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I].conjugate() * B[J][L];
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    }
  } else {
    if (CONJB) {
      // Form  C := alpha*A**T*B**H + beta*C

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I] * B[J][L].conjugate();
          }
          if (BETA == Complex.zero) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    } else {
      // Form  C := alpha*A**T*B**T + beta*C

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP += A[L][I] * B[J][L];
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
