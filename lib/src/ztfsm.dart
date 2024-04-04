import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void ztfsm(
  final String TRANSR,
  final String SIDE,
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int M,
  final int N,
  final Complex ALPHA,
  final Array<Complex> A_,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(offset: zeroIndexedArrayOffset);
  final B = B_.having(ld: LDB, offset: zeroIndexedMatrixOffset);
  int M1 = 0, M2 = 0, N1 = 0, N2 = 0, K = 0, INFO, I, J;

  // Test the input parameters.

  INFO = 0;
  final NORMALTRANSR = lsame(TRANSR, 'N');
  final LSIDE = lsame(SIDE, 'L');
  final LOWER = lsame(UPLO, 'L');
  final NOTRANS = lsame(TRANS, 'N');
  if (!NORMALTRANSR && !lsame(TRANSR, 'C')) {
    INFO = -1;
  } else if (!LSIDE && !lsame(SIDE, 'R')) {
    INFO = -2;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO = -3;
  } else if (!NOTRANS && !lsame(TRANS, 'C')) {
    INFO = -4;
  } else if (!lsame(DIAG, 'N') && !lsame(DIAG, 'U')) {
    INFO = -5;
  } else if (M < 0) {
    INFO = -6;
  } else if (N < 0) {
    INFO = -7;
  } else if (LDB < max(1, M)) {
    INFO = -11;
  }
  if (INFO != 0) {
    xerbla('ZTFSM', -INFO);
    return;
  }

  // Quick return when ( (N == 0) || (M == 0) )

  if ((M == 0) || (N == 0)) return;

  // Quick return when ALPHA == (0D+0,0D+0)

  if (ALPHA == Complex.zero) {
    for (J = 0; J <= N - 1; J++) {
      for (I = 0; I <= M - 1; I++) {
        B[I][J] = Complex.zero;
      }
    }
    return;
  }

  if (LSIDE) {
    // SIDE = 'L'

    // A is M-by-M.
    // If M is odd, set NISODD = true , and M1 and M2.
    // If M is even, NISODD = false , and M.
    bool MISODD;
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
              ztrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A.asMatrix(), M, B, LDB);
            } else {
              ztrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A(0).asMatrix(), M, B,
                  LDB);
              zgemm('N', 'N', M2, N, M1, -Complex.one, A(M1).asMatrix(), M, B,
                  LDB, ALPHA, B(M1, 0), LDB);
              ztrsm('L', 'U', 'C', DIAG, M2, N, Complex.one, A(M).asMatrix(), M,
                  B(M1, 0), LDB);
            }
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
            // TRANS = 'C'

            if (M == 1) {
              ztrsm('L', 'L', 'C', DIAG, M1, N, ALPHA, A(0).asMatrix(), M, B,
                  LDB);
            } else {
              ztrsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A(M).asMatrix(), M,
                  B(M1, 0), LDB);
              zgemm('C', 'N', M1, N, M2, -Complex.one, A(M1).asMatrix(), M,
                  B(M1, 0), LDB, ALPHA, B, LDB);
              ztrsm('L', 'L', 'C', DIAG, M1, N, Complex.one, A(0).asMatrix(), M,
                  B, LDB);
            }
          }
        } else {
          // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'N'

            ztrsm(
                'L', 'L', 'N', DIAG, M1, N, ALPHA, A(M2).asMatrix(), M, B, LDB);
            zgemm('C', 'N', M2, N, M1, -Complex.one, A(0).asMatrix(), M, B, LDB,
                ALPHA, B(M1, 0), LDB);
            ztrsm('L', 'U', 'C', DIAG, M2, N, Complex.one, A(M1).asMatrix(), M,
                B(M1, 0), LDB);
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'C'

            ztrsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A(M1).asMatrix(), M,
                B(M1, 0), LDB);
            zgemm('N', 'N', M1, N, M2, -Complex.one, A(0).asMatrix(), M,
                B(M1, 0), LDB, ALPHA, B, LDB);
            ztrsm('L', 'L', 'C', DIAG, M1, N, Complex.one, A(M2).asMatrix(), M,
                B, LDB);
          }
        }
      } else {
        // SIDE = 'L', N is odd, and TRANSR = 'C'

        if (LOWER) {
          // SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
            // TRANS = 'N'

            if (M == 1) {
              ztrsm('L', 'U', 'C', DIAG, M1, N, ALPHA, A(0).asMatrix(), M1, B,
                  LDB);
            } else {
              ztrsm('L', 'U', 'C', DIAG, M1, N, ALPHA, A(0).asMatrix(), M1, B,
                  LDB);
              zgemm('C', 'N', M2, N, M1, -Complex.one, A(M1 * M1).asMatrix(),
                  M1, B, LDB, ALPHA, B(M1, 0), LDB);
              ztrsm('L', 'L', 'N', DIAG, M2, N, Complex.one, A(1).asMatrix(),
                  M1, B(M1, 0), LDB);
            }
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
            // TRANS = 'C'

            if (M == 1) {
              ztrsm('L', 'U', 'N', DIAG, M1, N, ALPHA, A(0).asMatrix(), M1, B,
                  LDB);
            } else {
              ztrsm('L', 'L', 'C', DIAG, M2, N, ALPHA, A(1).asMatrix(), M1,
                  B(M1, 0), LDB);
              zgemm('N', 'N', M1, N, M2, -Complex.one, A(M1 * M1).asMatrix(),
                  M1, B(M1, 0), LDB, ALPHA, B, LDB);
              ztrsm('L', 'U', 'N', DIAG, M1, N, Complex.one, A(0).asMatrix(),
                  M1, B, LDB);
            }
          }
        } else {
          // SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
            // TRANS = 'N'

            ztrsm('L', 'U', 'C', DIAG, M1, N, ALPHA, A(M2 * M2).asMatrix(), M2,
                B, LDB);
            zgemm('N', 'N', M2, N, M1, -Complex.one, A(0).asMatrix(), M2, B,
                LDB, ALPHA, B(M1, 0), LDB);
            ztrsm('L', 'L', 'N', DIAG, M2, N, Complex.one,
                A(M1 * M2).asMatrix(), M2, B(M1, 0), LDB);
          } else {
            // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
            // TRANS = 'C'

            ztrsm('L', 'L', 'C', DIAG, M2, N, ALPHA, A(M1 * M2).asMatrix(), M2,
                B(M1, 0), LDB);
            zgemm('C', 'N', M1, N, M2, -Complex.one, A(0).asMatrix(), M2,
                B(M1, 0), LDB, ALPHA, B, LDB);
            ztrsm('L', 'U', 'N', DIAG, M1, N, Complex.one,
                A(M2 * M2).asMatrix(), M2, B, LDB);
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

            ztrsm('L', 'L', 'N', DIAG, K, N, ALPHA, A(1).asMatrix(), M + 1, B,
                LDB);
            zgemm('N', 'N', K, N, K, -Complex.one, A(K + 1).asMatrix(), M + 1,
                B, LDB, ALPHA, B(K, 0), LDB);
            ztrsm('L', 'U', 'C', DIAG, K, N, Complex.one, A(0).asMatrix(),
                M + 1, B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
            // and TRANS = 'C'

            ztrsm('L', 'U', 'N', DIAG, K, N, ALPHA, A(0).asMatrix(), M + 1,
                B(K, 0), LDB);
            zgemm('C', 'N', K, N, K, -Complex.one, A(K + 1).asMatrix(), M + 1,
                B(K, 0), LDB, ALPHA, B, LDB);
            ztrsm('L', 'L', 'C', DIAG, K, N, Complex.one, A(1).asMatrix(),
                M + 1, B, LDB);
          }
        } else {
          // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'N'

            ztrsm('L', 'L', 'N', DIAG, K, N, ALPHA, A(K + 1).asMatrix(), M + 1,
                B, LDB);
            zgemm('C', 'N', K, N, K, -Complex.one, A(0).asMatrix(), M + 1, B,
                LDB, ALPHA, B(K, 0), LDB);
            ztrsm('L', 'U', 'C', DIAG, K, N, Complex.one, A(K).asMatrix(),
                M + 1, B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'C'
            ztrsm('L', 'U', 'N', DIAG, K, N, ALPHA, A(K).asMatrix(), M + 1,
                B(K, 0), LDB);
            zgemm('N', 'N', K, N, K, -Complex.one, A(0).asMatrix(), M + 1,
                B(K, 0), LDB, ALPHA, B, LDB);
            ztrsm('L', 'L', 'C', DIAG, K, N, Complex.one, A(K + 1).asMatrix(),
                M + 1, B, LDB);
          }
        }
      } else {
        // SIDE = 'L', N is even, and TRANSR = 'C'

        if (LOWER) {
          // SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
            // and TRANS = 'N'

            ztrsm('L', 'U', 'C', DIAG, K, N, ALPHA, A(K).asMatrix(), K, B, LDB);
            zgemm('C', 'N', K, N, K, -Complex.one, A(K * (K + 1)).asMatrix(), K,
                B, LDB, ALPHA, B(K, 0), LDB);
            ztrsm('L', 'L', 'N', DIAG, K, N, Complex.one, A(0).asMatrix(), K,
                B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
            // and TRANS = 'C'

            ztrsm('L', 'L', 'C', DIAG, K, N, ALPHA, A(0).asMatrix(), K, B(K, 0),
                LDB);
            zgemm('N', 'N', K, N, K, -Complex.one, A(K * (K + 1)).asMatrix(), K,
                B(K, 0), LDB, ALPHA, B, LDB);
            ztrsm('L', 'U', 'N', DIAG, K, N, Complex.one, A(K).asMatrix(), K, B,
                LDB);
          }
        } else {
          // SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'U'

          if (!NOTRANS) {
            // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
            // and TRANS = 'N'

            ztrsm('L', 'U', 'C', DIAG, K, N, ALPHA, A(K * (K + 1)).asMatrix(),
                K, B, LDB);
            zgemm('N', 'N', K, N, K, -Complex.one, A(0).asMatrix(), K, B, LDB,
                ALPHA, B(K, 0), LDB);
            ztrsm('L', 'L', 'N', DIAG, K, N, Complex.one, A(K * K).asMatrix(),
                K, B(K, 0), LDB);
          } else {
            // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
            // and TRANS = 'C'

            ztrsm('L', 'L', 'C', DIAG, K, N, ALPHA, A(K * K).asMatrix(), K,
                B(K, 0), LDB);
            zgemm('C', 'N', K, N, K, -Complex.one, A(0).asMatrix(), K, B(K, 0),
                LDB, ALPHA, B, LDB);
            ztrsm('L', 'U', 'N', DIAG, K, N, Complex.one,
                A(K * (K + 1)).asMatrix(), K, B, LDB);
          }
        }
      }
    }
  } else {
    // SIDE = 'R'

    // A is N-by-N.
    // If N is odd, set NISODD = true , and N1 and N2.
    // If N is even, NISODD = false , and K.
    bool NISODD;
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

            ztrsm('R', 'U', 'C', DIAG, M, N2, ALPHA, A(N).asMatrix(), N,
                B(0, N1), LDB);
            zgemm('N', 'N', M, N1, N2, -Complex.one, B(0, N1), LDB,
                A(N1).asMatrix(), N, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'L', 'N', DIAG, M, N1, Complex.one, A(0).asMatrix(), N,
                B(0, 0), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
            // TRANS = 'C'

            ztrsm('R', 'L', 'C', DIAG, M, N1, ALPHA, A(0).asMatrix(), N,
                B(0, 0), LDB);
            zgemm('N', 'C', M, N2, N1, -Complex.one, B(0, 0), LDB,
                A(N1).asMatrix(), N, ALPHA, B(0, N1), LDB);
            ztrsm('R', 'U', 'N', DIAG, M, N2, Complex.one, A(N).asMatrix(), N,
                B(0, N1), LDB);
          }
        } else {
          // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'N'

            ztrsm('R', 'L', 'C', DIAG, M, N1, ALPHA, A(N2).asMatrix(), N,
                B(0, 0), LDB);
            zgemm('N', 'N', M, N2, N1, -Complex.one, B(0, 0), LDB,
                A(0).asMatrix(), N, ALPHA, B(0, N1), LDB);
            ztrsm('R', 'U', 'N', DIAG, M, N2, Complex.one, A(N1).asMatrix(), N,
                B(0, N1), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
            // TRANS = 'C'

            ztrsm('R', 'U', 'C', DIAG, M, N2, ALPHA, A(N1).asMatrix(), N,
                B(0, N1), LDB);
            zgemm('N', 'C', M, N1, N2, -Complex.one, B(0, N1), LDB,
                A(0).asMatrix(), N, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'L', 'N', DIAG, M, N1, Complex.one, A(N2).asMatrix(), N,
                B(0, 0), LDB);
          }
        }
      } else {
        // SIDE = 'R', N is odd, and TRANSR = 'C'

        if (LOWER) {
          // SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
            // TRANS = 'N'

            ztrsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A(1).asMatrix(), N1,
                B(0, N1), LDB);
            zgemm('N', 'C', M, N1, N2, -Complex.one, B(0, N1), LDB,
                A(N1 * N1).asMatrix(), N1, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'U', 'C', DIAG, M, N1, Complex.one, A(0).asMatrix(), N1,
                B(0, 0), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
            // TRANS = 'C'

            ztrsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A(0).asMatrix(), N1,
                B(0, 0), LDB);
            zgemm('N', 'N', M, N2, N1, -Complex.one, B(0, 0), LDB,
                A(N1 * N1).asMatrix(), N1, ALPHA, B(0, N1), LDB);
            ztrsm('R', 'L', 'C', DIAG, M, N2, Complex.one, A(1).asMatrix(), N1,
                B(0, N1), LDB);
          }
        } else {
          // SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
            // TRANS = 'N'

            ztrsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A(N2 * N2).asMatrix(), N2,
                B(0, 0), LDB);
            zgemm('N', 'C', M, N2, N1, -Complex.one, B(0, 0), LDB,
                A(0).asMatrix(), N2, ALPHA, B(0, N1), LDB);
            ztrsm('R', 'L', 'C', DIAG, M, N2, Complex.one,
                A(N1 * N2).asMatrix(), N2, B(0, N1), LDB);
          } else {
            // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
            // TRANS = 'C'

            ztrsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A(N1 * N2).asMatrix(), N2,
                B(0, N1), LDB);
            zgemm('N', 'N', M, N1, N2, -Complex.one, B(0, N1), LDB,
                A(0).asMatrix(), N2, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'U', 'C', DIAG, M, N1, Complex.one,
                A(N2 * N2).asMatrix(), N2, B(0, 0), LDB);
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

            ztrsm('R', 'U', 'C', DIAG, M, K, ALPHA, A(0).asMatrix(), N + 1,
                B(0, K), LDB);
            zgemm('N', 'N', M, K, K, -Complex.one, B(0, K), LDB,
                A(K + 1).asMatrix(), N + 1, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'L', 'N', DIAG, M, K, Complex.one, A(1).asMatrix(),
                N + 1, B(0, 0), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
            // and TRANS = 'C'

            ztrsm('R', 'L', 'C', DIAG, M, K, ALPHA, A(1).asMatrix(), N + 1,
                B(0, 0), LDB);
            zgemm('N', 'C', M, K, K, -Complex.one, B(0, 0), LDB,
                A(K + 1).asMatrix(), N + 1, ALPHA, B(0, K), LDB);
            ztrsm('R', 'U', 'N', DIAG, M, K, Complex.one, A(0).asMatrix(),
                N + 1, B(0, K), LDB);
          }
        } else {
          // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'N'

            ztrsm('R', 'L', 'C', DIAG, M, K, ALPHA, A(K + 1).asMatrix(), N + 1,
                B(0, 0), LDB);
            zgemm('N', 'N', M, K, K, -Complex.one, B(0, 0), LDB,
                A(0).asMatrix(), N + 1, ALPHA, B(0, K), LDB);
            ztrsm('R', 'U', 'N', DIAG, M, K, Complex.one, A(K).asMatrix(),
                N + 1, B(0, K), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
            // and TRANS = 'C'

            ztrsm('R', 'U', 'C', DIAG, M, K, ALPHA, A(K).asMatrix(), N + 1,
                B(0, K), LDB);
            zgemm('N', 'C', M, K, K, -Complex.one, B(0, K), LDB,
                A(0).asMatrix(), N + 1, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'L', 'N', DIAG, M, K, Complex.one, A(K + 1).asMatrix(),
                N + 1, B(0, 0), LDB);
          }
        }
      } else {
        // SIDE = 'R', N is even, and TRANSR = 'C'

        if (LOWER) {
          // SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'L'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
            // and TRANS = 'N'

            ztrsm('R', 'L', 'N', DIAG, M, K, ALPHA, A(0).asMatrix(), K, B(0, K),
                LDB);
            zgemm('N', 'C', M, K, K, -Complex.one, B(0, K), LDB,
                A((K + 1) * K).asMatrix(), K, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'U', 'C', DIAG, M, K, Complex.one, A(K).asMatrix(), K,
                B(0, 0), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
            // and TRANS = 'C'

            ztrsm('R', 'U', 'N', DIAG, M, K, ALPHA, A(K).asMatrix(), K, B(0, 0),
                LDB);
            zgemm('N', 'N', M, K, K, -Complex.one, B(0, 0), LDB,
                A((K + 1) * K).asMatrix(), K, ALPHA, B(0, K), LDB);
            ztrsm('R', 'L', 'C', DIAG, M, K, Complex.one, A(0).asMatrix(), K,
                B(0, K), LDB);
          }
        } else {
          // SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'U'

          if (NOTRANS) {
            // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
            // and TRANS = 'N'

            ztrsm('R', 'U', 'N', DIAG, M, K, ALPHA, A((K + 1) * K).asMatrix(),
                K, B(0, 0), LDB);
            zgemm('N', 'C', M, K, K, -Complex.one, B(0, 0), LDB,
                A(0).asMatrix(), K, ALPHA, B(0, K), LDB);
            ztrsm('R', 'L', 'C', DIAG, M, K, Complex.one, A(K * K).asMatrix(),
                K, B(0, K), LDB);
          } else {
            // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
            // and TRANS = 'C'

            ztrsm('R', 'L', 'N', DIAG, M, K, ALPHA, A(K * K).asMatrix(), K,
                B(0, K), LDB);
            zgemm('N', 'N', M, K, K, -Complex.one, B(0, K), LDB,
                A(0).asMatrix(), K, ALPHA, B(0, 0), LDB);
            ztrsm('R', 'U', 'C', DIAG, M, K, Complex.one,
                A((K + 1) * K).asMatrix(), K, B(0, 0), LDB);
          }
        }
      }
    }
  }
}
