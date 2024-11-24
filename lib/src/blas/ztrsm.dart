// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void ztrsm(
  final String SIDE,
  final String UPLO,
  final String TRANSA,
  final String DIAG,
  final int M,
  final int N,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  Complex TEMP;
  int I, INFO, J, K, NROWA;
  bool LSIDE, NOCONJ, NOUNIT, UPPER;

  // Test the input parameters.

  LSIDE = lsame(SIDE, 'L');
  if (LSIDE) {
    NROWA = M;
  } else {
    NROWA = N;
  }
  NOCONJ = lsame(TRANSA, 'T');
  NOUNIT = lsame(DIAG, 'N');
  UPPER = lsame(UPLO, 'U');

  INFO = 0;
  if (!LSIDE && !lsame(SIDE, 'R')) {
    INFO = 1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO = 2;
  } else if (!lsame(TRANSA, 'N') &&
      !lsame(TRANSA, 'T') &&
      !lsame(TRANSA, 'C')) {
    INFO = 3;
  } else if (!lsame(DIAG, 'U') && !lsame(DIAG, 'N')) {
    INFO = 4;
  } else if (M < 0) {
    INFO = 5;
  } else if (N < 0) {
    INFO = 6;
  } else if (LDA < max(1, NROWA)) {
    INFO = 9;
  } else if (LDB < max(1, M)) {
    INFO = 11;
  }
  if (INFO != 0) {
    xerbla('ZTRSM', INFO);
    return;
  }

  // Quick return if possible.

  if (M == 0 || N == 0) return;

  // And when  alpha == zero.

  if (ALPHA == Complex.zero) {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= M; I++) {
        B[I][J] = Complex.zero;
      }
    }
    return;
  }

  // Start the operations.

  if (LSIDE) {
    if (lsame(TRANSA, 'N')) {
      // Form  B := alpha*inv( A )*B.

      if (UPPER) {
        for (J = 1; J <= N; J++) {
          if (ALPHA != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][J] = ALPHA * B[I][J];
            }
          }
          for (K = M; K >= 1; K--) {
            if (B[K][J] != Complex.zero) {
              if (NOUNIT) B[K][J] /= A[K][K];
              for (I = 1; I <= K - 1; I++) {
                B[I][J] -= B[K][J] * A[I][K];
              }
            }
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          if (ALPHA != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][J] = ALPHA * B[I][J];
            }
          }
          for (K = 1; K <= M; K++) {
            if (B[K][J] != Complex.zero) {
              if (NOUNIT) B[K][J] /= A[K][K];
              for (I = K + 1; I <= M; I++) {
                B[I][J] -= B[K][J] * A[I][K];
              }
            }
          }
        }
      }
    } else {
      // Form  B := alpha*inv( A**T )*B
      // or    B := alpha*inv( A**H )*B.

      if (UPPER) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= M; I++) {
            TEMP = ALPHA * B[I][J];
            if (NOCONJ) {
              for (K = 1; K <= I - 1; K++) {
                TEMP -= A[K][I] * B[K][J];
              }
              if (NOUNIT) TEMP /= A[I][I];
            } else {
              for (K = 1; K <= I - 1; K++) {
                TEMP -= A[K][I].conjugate() * B[K][J];
              }
              if (NOUNIT) TEMP /= A[I][I].conjugate();
            }
            B[I][J] = TEMP;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = M; I >= 1; I--) {
            TEMP = ALPHA * B[I][J];
            if (NOCONJ) {
              for (K = I + 1; K <= M; K++) {
                TEMP -= A[K][I] * B[K][J];
              }
              if (NOUNIT) TEMP /= A[I][I];
            } else {
              for (K = I + 1; K <= M; K++) {
                TEMP -= A[K][I].conjugate() * B[K][J];
              }
              if (NOUNIT) TEMP /= A[I][I].conjugate();
            }
            B[I][J] = TEMP;
          }
        }
      }
    }
  } else {
    if (lsame(TRANSA, 'N')) {
      // Form  B := alpha*B*inv( A ).

      if (UPPER) {
        for (J = 1; J <= N; J++) {
          if (ALPHA != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][J] = ALPHA * B[I][J];
            }
          }
          for (K = 1; K <= J - 1; K++) {
            if (A[K][J] != Complex.zero) {
              for (I = 1; I <= M; I++) {
                B[I][J] -= A[K][J] * B[I][K];
              }
            }
          }
          if (NOUNIT) {
            TEMP = Complex.one / A[J][J];
            for (I = 1; I <= M; I++) {
              B[I][J] = TEMP * B[I][J];
            }
          }
        }
      } else {
        for (J = N; J >= 1; J--) {
          if (ALPHA != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][J] = ALPHA * B[I][J];
            }
          }
          for (K = J + 1; K <= N; K++) {
            if (A[K][J] != Complex.zero) {
              for (I = 1; I <= M; I++) {
                B[I][J] -= A[K][J] * B[I][K];
              }
            }
          }
          if (NOUNIT) {
            TEMP = Complex.one / A[J][J];
            for (I = 1; I <= M; I++) {
              B[I][J] = TEMP * B[I][J];
            }
          }
        }
      }
    } else {
      // Form  B := alpha*B*inv( A**T )
      // or    B := alpha*B*inv( A**H ).

      if (UPPER) {
        for (K = N; K >= 1; K--) {
          if (NOUNIT) {
            if (NOCONJ) {
              TEMP = Complex.one / A[K][K];
            } else {
              TEMP = Complex.one / A[K][K].conjugate();
            }
            for (I = 1; I <= M; I++) {
              B[I][K] = TEMP * B[I][K];
            }
          }
          for (J = 1; J <= K - 1; J++) {
            if (A[J][K] != Complex.zero) {
              if (NOCONJ) {
                TEMP = A[J][K];
              } else {
                TEMP = A[J][K].conjugate();
              }
              for (I = 1; I <= M; I++) {
                B[I][J] -= TEMP * B[I][K];
              }
            }
          }
          if (ALPHA != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][K] = ALPHA * B[I][K];
            }
          }
        }
      } else {
        for (K = 1; K <= N; K++) {
          if (NOUNIT) {
            if (NOCONJ) {
              TEMP = Complex.one / A[K][K];
            } else {
              TEMP = Complex.one / A[K][K].conjugate();
            }
            for (I = 1; I <= M; I++) {
              B[I][K] = TEMP * B[I][K];
            }
          }
          for (J = K + 1; J <= N; J++) {
            if (A[J][K] != Complex.zero) {
              if (NOCONJ) {
                TEMP = A[J][K];
              } else {
                TEMP = A[J][K].conjugate();
              }
              for (I = 1; I <= M; I++) {
                B[I][J] -= TEMP * B[I][K];
              }
            }
          }
          if (ALPHA != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][K] = ALPHA * B[I][K];
            }
          }
        }
      }
    }
  }
}
