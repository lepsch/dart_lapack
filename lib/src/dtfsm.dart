// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtfsm(
  final String TRANSR,
  final String SIDE,
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int M,
  final int N,
  final double ALPHA,
  final Array<double> A_,
  final Matrix<double> B_,
  final int LDB,
) {
  final A = A_.having(offset: zeroIndexedArrayOffset);
  final B = B_.having(ld: LDB, offset: zeroIndexedMatrixOffset);
  const ONE = 1.0, ZERO = 0.0;
  bool LOWER, LSIDE, MISODD, NISODD, NORMALTRANSR, NOTRANS;
  int M1 = 0, M2 = 0, N1 = 0, N2 = 0, K = 0, I, J;
  final INFO = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LSIDE = lsame(SIDE, 'L');
  LOWER = lsame(UPLO, 'L');
  NOTRANS = lsame(TRANS, 'N');
  if (!NORMALTRANSR && !lsame(TRANSR, 'T')) {
    INFO.value = -1;
  } else if (!LSIDE && !lsame(SIDE, 'R')) {
    INFO.value = -2;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -3;
  } else if (!NOTRANS && !lsame(TRANS, 'T')) {
    INFO.value = -4;
  } else if (!lsame(DIAG, 'N') && !lsame(DIAG, 'U')) {
    INFO.value = -5;
  } else if (M < 0) {
    INFO.value = -6;
  } else if (N < 0) {
    INFO.value = -7;
  } else if (LDB < max(1, M)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DTFSM', -INFO.value);
    return;
  }

  // Quick return when ( (N == 0) || (M == 0) )

  if ((M == 0) || (N == 0)) return;

  // Quick return when ALPHA == (0D+0)

  if (ALPHA == ZERO) {
    for (J = 0; J <= N - 1; J++) {
      for (I = 0; I <= M - 1; I++) {
        B[I][J] = ZERO;
      }
    }
    return;
  }

  if (LSIDE) {
    // SIDE = 'L'

    // A is M-by-M.
    // If M is odd, set NISODD = true , and M1 and M2.
    // If M is even, NISODD = false , and M.

    if ((M % 2) == 0) {
      MISODD = false;
      K = M ~/ 2;
    } else {
      MISODD = true;
      if (LOWER) {
        M2 = M ~/ 2;
        M1 = M - M2;
      } else {
        M1 = M ~/ 2;
        M2 = M - M1;
      }
    }

    if (MISODD) {
      // SIDE = 'L' and N is odd

      if (NORMALTRANSR) {
        // SIDE = 'L', N is odd, and TRANSR = 'N'

        if (LOWER) {
          // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
            // TRANS = 'N'

            if (M == 1) {
              dtrsm(
                  'L', 'L', 'N', DIAG, M1, N, ALPHA, A.asMatrix(M), M, B, LDB);
            } else {
              dtrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A(0).asMatrix(M), M, B,
                  LDB);
              dgemm('N', 'N', M2, N, M1, -ONE, A(M1).asMatrix(M), M, B, LDB,
                  ALPHA, B(M1, 0), LDB);
              dtrsm('L', 'U', 'T', DIAG, M2, N, ONE, A(M).asMatrix(M), M,
                  B(M1, 0), LDB);
            }
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
            // TRANS = 'T'

            if (M == 1) {
              dtrsm('L', 'L', 'T', DIAG, M1, N, ALPHA, A(0).asMatrix(M), M, B,
                  LDB);
            } else {
              dtrsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A(M).asMatrix(M), M,
                  B(M1, 0), LDB);
              dgemm('T', 'N', M1, N, M2, -ONE, A(M1).asMatrix(M), M, B(M1, 0),
                  LDB, ALPHA, B, LDB);
              dtrsm(
                  'L', 'L', 'T', DIAG, M1, N, ONE, A(0).asMatrix(M), M, B, LDB);
            }
          }
        } else {
          // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'N'

            dtrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A(M2).asMatrix(M), M, B,
                LDB);
            dgemm('T', 'N', M2, N, M1, -ONE, A(0).asMatrix(M), M, B, LDB, ALPHA,
                B(M1, 0), LDB);
            dtrsm('L', 'U', 'T', DIAG, M2, N, ONE, A(M1).asMatrix(M), M,
                B(M1, 0), LDB);
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'T'

            dtrsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A(M1).asMatrix(M), M,
                B(M1, 0), LDB);
            dgemm('N', 'N', M1, N, M2, -ONE, A(0).asMatrix(M), M, B(M1, 0), LDB,
                ALPHA, B, LDB);
            dtrsm(
                'L', 'L', 'T', DIAG, M1, N, ONE, A(M2).asMatrix(M), M, B, LDB);
          }
        }
      } else {
        // SIDE = 'L', N is odd, and TRANSR = 'T'

        if (LOWER) {
          // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
            // TRANS = 'N'

            if (M == 1) {
              dtrsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A(0).asMatrix(M1), M1, B,
                  LDB);
            } else {
              dtrsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A(0).asMatrix(M1), M1, B,
                  LDB);
              dgemm('T', 'N', M2, N, M1, -ONE, A(M1 * M1).asMatrix(M1), M1, B,
                  LDB, ALPHA, B(M1, 0), LDB);
              dtrsm('L', 'L', 'N', DIAG, M2, N, ONE, A(1).asMatrix(M1), M1,
                  B(M1, 0), LDB);
            }
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
            // TRANS = 'T'

            if (M == 1) {
              dtrsm('L', 'U', 'N', DIAG, M1, N, ALPHA, A(0).asMatrix(M1), M1, B,
                  LDB);
            } else {
              dtrsm('L', 'L', 'T', DIAG, M2, N, ALPHA, A(1).asMatrix(M1), M1,
                  B(M1, 0), LDB);
              dgemm('N', 'N', M1, N, M2, -ONE, A(M1 * M1).asMatrix(M1), M1,
                  B(M1, 0), LDB, ALPHA, B, LDB);
              dtrsm('L', 'U', 'N', DIAG, M1, N, ONE, A(0).asMatrix(M1), M1, B,
                  LDB);
            }
          }
        } else {
          // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
            // TRANS = 'N'

            dtrsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A(M2 * M2).asMatrix(M2),
                M2, B, LDB);
            dgemm('N', 'N', M2, N, M1, -ONE, A(0).asMatrix(M2), M2, B, LDB,
                ALPHA, B(M1, 0), LDB);
            dtrsm('L', 'L', 'N', DIAG, M2, N, ONE, A(M1 * M2).asMatrix(M2), M2,
                B(M1, 0), LDB);
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
            // TRANS = 'T'

            dtrsm('L', 'L', 'T', DIAG, M2, N, ALPHA, A(M1 * M2).asMatrix(M2),
                M2, B(M1, 0), LDB);
            dgemm('T', 'N', M1, N, M2, -ONE, A(0).asMatrix(M2), M2, B(M1, 0),
                LDB, ALPHA, B, LDB);
            dtrsm('L', 'U', 'N', DIAG, M1, N, ONE, A(M2 * M2).asMatrix(M2), M2,
                B, LDB);
          }
        }
      }
    } else {
      // SIDE = 'L' and N is even

      if (NORMALTRANSR) {
        // SIDE = 'L', N is even, and TRANSR = 'N'

        if (LOWER) {
          // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
            // and TRANS = 'N'

            dtrsm('L', 'L', 'N', DIAG, K, N, ALPHA, A(1).asMatrix(M + 1), M + 1,
                B, LDB);
            dgemm('N', 'N', K, N, K, -ONE, A(K + 1).asMatrix(M + 1), M + 1, B,
                LDB, ALPHA, B(K, 0), LDB);
            dtrsm('L', 'U', 'T', DIAG, K, N, ONE, A(0).asMatrix(M + 1), M + 1,
                B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
            // and TRANS = 'T'

            dtrsm('L', 'U', 'N', DIAG, K, N, ALPHA, A(0).asMatrix(M + 1), M + 1,
                B(K, 0), LDB);
            dgemm('T', 'N', K, N, K, -ONE, A(K + 1).asMatrix(M + 1), M + 1,
                B(K, 0), LDB, ALPHA, B, LDB);
            dtrsm('L', 'L', 'T', DIAG, K, N, ONE, A(1).asMatrix(M + 1), M + 1,
                B, LDB);
          }
        } else {
          // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'N'

            dtrsm('L', 'L', 'N', DIAG, K, N, ALPHA, A(K + 1).asMatrix(M + 1),
                M + 1, B, LDB);
            dgemm('T', 'N', K, N, K, -ONE, A(0).asMatrix(M + 1), M + 1, B, LDB,
                ALPHA, B(K, 0), LDB);
            dtrsm('L', 'U', 'T', DIAG, K, N, ONE, A(K).asMatrix(M + 1), M + 1,
                B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'T'
            dtrsm('L', 'U', 'N', DIAG, K, N, ALPHA, A(K).asMatrix(M + 1), M + 1,
                B(K, 0), LDB);
            dgemm('N', 'N', K, N, K, -ONE, A(0).asMatrix(M + 1), M + 1, B(K, 0),
                LDB, ALPHA, B, LDB);
            dtrsm('L', 'L', 'T', DIAG, K, N, ONE, A(K + 1).asMatrix(M + 1),
                M + 1, B, LDB);
          }
        }
      } else {
        // SIDE = 'L', N is even, and TRANSR = 'T'

        if (LOWER) {
          // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
            // and TRANS = 'N'

            dtrsm(
                'L', 'U', 'T', DIAG, K, N, ALPHA, A(K).asMatrix(K), K, B, LDB);
            dgemm('T', 'N', K, N, K, -ONE, A(K * (K + 1)).asMatrix(K), K, B,
                LDB, ALPHA, B(K, 0), LDB);
            dtrsm('L', 'L', 'N', DIAG, K, N, ONE, A(0).asMatrix(K), K, B(K, 0),
                LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
            // and TRANS = 'T'

            dtrsm('L', 'L', 'T', DIAG, K, N, ALPHA, A(0).asMatrix(K), K,
                B(K, 0), LDB);
            dgemm('N', 'N', K, N, K, -ONE, A(K * (K + 1)).asMatrix(K), K,
                B(K, 0), LDB, ALPHA, B, LDB);
            dtrsm('L', 'U', 'N', DIAG, K, N, ONE, A(K).asMatrix(K), K, B, LDB);
          }
        } else {
          // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
            // and TRANS = 'N'

            dtrsm('L', 'U', 'T', DIAG, K, N, ALPHA, A(K * (K + 1)).asMatrix(K),
                K, B, LDB);
            dgemm('N', 'N', K, N, K, -ONE, A(0).asMatrix(K), K, B, LDB, ALPHA,
                B(K, 0), LDB);
            dtrsm('L', 'L', 'N', DIAG, K, N, ONE, A(K * K).asMatrix(K), K,
                B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
            // and TRANS = 'T'

            dtrsm('L', 'L', 'T', DIAG, K, N, ALPHA, A(K * K).asMatrix(K), K,
                B(K, 0), LDB);
            dgemm('T', 'N', K, N, K, -ONE, A(0).asMatrix(K), K, B(K, 0), LDB,
                ALPHA, B, LDB);
            dtrsm('L', 'U', 'N', DIAG, K, N, ONE, A(K * (K + 1)).asMatrix(K), K,
                B, LDB);
          }
        }
      }
    }
  } else {
    // SIDE = 'R'

    // A is N-by-N.
    // If N is odd, set NISODD = true , and N1 and N2.
    // If N is even, NISODD = false , and K.

    if ((N % 2) == 0) {
      NISODD = false;
      K = N ~/ 2;
    } else {
      NISODD = true;
      if (LOWER) {
        N2 = N ~/ 2;
        N1 = N - N2;
      } else {
        N1 = N ~/ 2;
        N2 = N - N1;
      }
    }

    if (NISODD) {
      // SIDE = 'R' and N is odd

      if (NORMALTRANSR) {
        // SIDE = 'R', N is odd, and TRANSR = 'N'

        if (LOWER) {
          // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
            // TRANS = 'N'

            dtrsm('R', 'U', 'T', DIAG, M, N2, ALPHA, A(N).asMatrix(N), N,
                B(0, N1), LDB);
            dgemm('N', 'N', M, N1, N2, -ONE, B(0, N1), LDB, A(N1).asMatrix(N),
                N, ALPHA, B(0, 0), LDB);
            dtrsm('R', 'L', 'N', DIAG, M, N1, ONE, A(0).asMatrix(N), N, B(0, 0),
                LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
            // TRANS = 'T'

            dtrsm('R', 'L', 'T', DIAG, M, N1, ALPHA, A(0).asMatrix(N), N,
                B(0, 0), LDB);
            dgemm('N', 'T', M, N2, N1, -ONE, B(0, 0), LDB, A(N1).asMatrix(N), N,
                ALPHA, B(0, N1), LDB);
            dtrsm('R', 'U', 'N', DIAG, M, N2, ONE, A(N).asMatrix(N), N,
                B(0, N1), LDB);
          }
        } else {
          // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'N'

            dtrsm('R', 'L', 'T', DIAG, M, N1, ALPHA, A(N2).asMatrix(N), N,
                B(0, 0), LDB);
            dgemm('N', 'N', M, N2, N1, -ONE, B(0, 0), LDB, A(0).asMatrix(N), N,
                ALPHA, B(0, N1), LDB);
            dtrsm('R', 'U', 'N', DIAG, M, N2, ONE, A(N1).asMatrix(N), N,
                B(0, N1), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'T'

            dtrsm('R', 'U', 'T', DIAG, M, N2, ALPHA, A(N1).asMatrix(N), N,
                B(0, N1), LDB);
            dgemm('N', 'T', M, N1, N2, -ONE, B(0, N1), LDB, A(0).asMatrix(N), N,
                ALPHA, B(0, 0), LDB);
            dtrsm('R', 'L', 'N', DIAG, M, N1, ONE, A(N2).asMatrix(N), N,
                B(0, 0), LDB);
          }
        }
      } else {
        // SIDE = 'R', N is odd, and TRANSR = 'T'

        if (LOWER) {
          // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
            // TRANS = 'N'

            dtrsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A(1).asMatrix(N1), N1,
                B(0, N1), LDB);
            dgemm('N', 'T', M, N1, N2, -ONE, B(0, N1), LDB,
                A(N1 * N1).asMatrix(N1), N1, ALPHA, B(0, 0), LDB);
            dtrsm('R', 'U', 'T', DIAG, M, N1, ONE, A(0).asMatrix(N1), N1,
                B(0, 0), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
            // TRANS = 'T'

            dtrsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A(0).asMatrix(N1), N1,
                B(0, 0), LDB);
            dgemm('N', 'N', M, N2, N1, -ONE, B(0, 0), LDB,
                A(N1 * N1).asMatrix(N1), N1, ALPHA, B(0, N1), LDB);
            dtrsm('R', 'L', 'T', DIAG, M, N2, ONE, A(1).asMatrix(N1), N1,
                B(0, N1), LDB);
          }
        } else {
          // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
            // TRANS = 'N'

            dtrsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A(N2 * N2).asMatrix(N2),
                N2, B(0, 0), LDB);
            dgemm('N', 'T', M, N2, N1, -ONE, B(0, 0), LDB, A(0).asMatrix(N2),
                N2, ALPHA, B(0, N1), LDB);
            dtrsm('R', 'L', 'T', DIAG, M, N2, ONE, A(N1 * N2).asMatrix(N2), N2,
                B(0, N1), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
            // TRANS = 'T'

            dtrsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A(N1 * N2).asMatrix(N2),
                N2, B(0, N1), LDB);
            dgemm('N', 'N', M, N1, N2, -ONE, B(0, N1), LDB, A(0).asMatrix(N2),
                N2, ALPHA, B(0, 0), LDB);
            dtrsm('R', 'U', 'T', DIAG, M, N1, ONE, A(N2 * N2).asMatrix(N2), N2,
                B(0, 0), LDB);
          }
        }
      }
    } else {
      // SIDE = 'R' and N is even

      if (NORMALTRANSR) {
        // SIDE = 'R', N is even, and TRANSR = 'N'

        if (LOWER) {
          // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
            // and TRANS = 'N'

            dtrsm('R', 'U', 'T', DIAG, M, K, ALPHA, A(0).asMatrix(N + 1), N + 1,
                B(0, K), LDB);
            dgemm('N', 'N', M, K, K, -ONE, B(0, K), LDB,
                A(K + 1).asMatrix(N + 1), N + 1, ALPHA, B(0, 0), LDB);
            dtrsm('R', 'L', 'N', DIAG, M, K, ONE, A(1).asMatrix(N + 1), N + 1,
                B(0, 0), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
            // and TRANS = 'T'

            dtrsm('R', 'L', 'T', DIAG, M, K, ALPHA, A(1).asMatrix(N + 1), N + 1,
                B(0, 0), LDB);
            dgemm('N', 'T', M, K, K, -ONE, B(0, 0), LDB,
                A(K + 1).asMatrix(N + 1), N + 1, ALPHA, B(0, K), LDB);
            dtrsm('R', 'U', 'N', DIAG, M, K, ONE, A(0).asMatrix(N + 1), N + 1,
                B(0, K), LDB);
          }
        } else {
          // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'N'

            dtrsm('R', 'L', 'T', DIAG, M, K, ALPHA, A(K + 1).asMatrix(N + 1),
                N + 1, B(0, 0), LDB);
            dgemm('N', 'N', M, K, K, -ONE, B(0, 0), LDB, A(0).asMatrix(N + 1),
                N + 1, ALPHA, B(0, K), LDB);
            dtrsm('R', 'U', 'N', DIAG, M, K, ONE, A(K).asMatrix(N + 1), N + 1,
                B(0, K), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'T'

            dtrsm('R', 'U', 'T', DIAG, M, K, ALPHA, A(K).asMatrix(N + 1), N + 1,
                B(0, K), LDB);
            dgemm('N', 'T', M, K, K, -ONE, B(0, K), LDB, A(0).asMatrix(N + 1),
                N + 1, ALPHA, B(0, 0), LDB);
            dtrsm('R', 'L', 'N', DIAG, M, K, ONE, A(K + 1).asMatrix(N + 1),
                N + 1, B(0, 0), LDB);
          }
        }
      } else {
        // SIDE = 'R', N is even, and TRANSR = 'T'

        if (LOWER) {
          // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
            // and TRANS = 'N'

            dtrsm('R', 'L', 'N', DIAG, M, K, ALPHA, A(0).asMatrix(K), K,
                B(0, K), LDB);
            dgemm('N', 'T', M, K, K, -ONE, B(0, K), LDB,
                A((K + 1) * K).asMatrix(K), K, ALPHA, B(0, 0), LDB);
            dtrsm('R', 'U', 'T', DIAG, M, K, ONE, A(K).asMatrix(K), K, B(0, 0),
                LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
            // and TRANS = 'T'

            dtrsm('R', 'U', 'N', DIAG, M, K, ALPHA, A(K).asMatrix(K), K,
                B(0, 0), LDB);
            dgemm('N', 'N', M, K, K, -ONE, B(0, 0), LDB,
                A((K + 1) * K).asMatrix(K), K, ALPHA, B(0, K), LDB);
            dtrsm('R', 'L', 'T', DIAG, M, K, ONE, A(0).asMatrix(K), K, B(0, K),
                LDB);
          }
        } else {
          // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
            // and TRANS = 'N'

            dtrsm('R', 'U', 'N', DIAG, M, K, ALPHA, A((K + 1) * K).asMatrix(K),
                K, B(0, 0), LDB);
            dgemm('N', 'T', M, K, K, -ONE, B(0, 0), LDB, A(0).asMatrix(K), K,
                ALPHA, B(0, K), LDB);
            dtrsm('R', 'L', 'T', DIAG, M, K, ONE, A(K * K).asMatrix(K), K,
                B(0, K), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
            // and TRANS = 'T'

            dtrsm('R', 'L', 'N', DIAG, M, K, ALPHA, A(K * K).asMatrix(K), K,
                B(0, K), LDB);
            dgemm('N', 'N', M, K, K, -ONE, B(0, K), LDB, A(0).asMatrix(K), K,
                ALPHA, B(0, 0), LDB);
            dtrsm('R', 'U', 'T', DIAG, M, K, ONE, A((K + 1) * K).asMatrix(K), K,
                B(0, 0), LDB);
          }
        }
      }
    }
  }
}
