import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void ztrmm(
  final String SIDE,
  final String UPLO,
  final String TRANSA,
  final String DIAG,
  final int M,
  final int N,
  final Complex ALPHA,
  final Matrix<Complex> A,
  final int LDA,
  final Matrix<Complex> B,
  final int LDB,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
  if ((!LSIDE) && (!lsame(SIDE, 'R'))) {
    INFO = 1;
  } else if ((!UPPER) && (!lsame(UPLO, 'L'))) {
    INFO = 2;
  } else if ((!lsame(TRANSA, 'N')) &&
      (!lsame(TRANSA, 'T')) &&
      (!lsame(TRANSA, 'C'))) {
    INFO = 3;
  } else if ((!lsame(DIAG, 'U')) && (!lsame(DIAG, 'N'))) {
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
    xerbla('ZTRMM ', INFO);
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
      // Form  B := alpha*A*B.

      if (UPPER) {
        for (J = 1; J <= N; J++) {
          for (K = 1; K <= M; K++) {
            if (B[K][J] != Complex.zero) {
              TEMP = ALPHA * B[K][J];
              for (I = 1; I <= K - 1; I++) {
                B[I][J] = B[I][J] + TEMP * A[I][K];
              }
              if (NOUNIT) TEMP = TEMP * A[K][K];
              B[K][J] = TEMP;
            }
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (K = M; K >= 1; K--) {
            if (B[K][J] != Complex.zero) {
              TEMP = ALPHA * B[K][J];
              B[K][J] = TEMP;
              if (NOUNIT) B[K][J] = B[K][J] * A[K][K];
              for (I = K + 1; I <= M; I++) {
                B[I][J] = B[I][J] + TEMP * A[I][K];
              }
            }
          }
        }
      }
    } else {
      // Form  B := alpha*A**T*B   or   B := alpha*A**H*B.

      if (UPPER) {
        for (J = 1; J <= N; J++) {
          for (I = M; I >= 1; I--) {
            TEMP = B[I][J];
            if (NOCONJ) {
              if (NOUNIT) TEMP = TEMP * A[I][I];
              for (K = 1; K <= I - 1; K++) {
                TEMP = TEMP + A[K][I] * B[K][J];
              }
            } else {
              if (NOUNIT) TEMP = TEMP * A[I][I].conjugate();
              for (K = 1; K <= I - 1; K++) {
                TEMP = TEMP + A[K][I].conjugate() * B[K][J];
              }
            }
            B[I][J] = ALPHA * TEMP;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= M; I++) {
            TEMP = B[I][J];
            if (NOCONJ) {
              if (NOUNIT) TEMP = TEMP * A[I][I];
              for (K = I + 1; K <= M; K++) {
                TEMP = TEMP + A[K][I] * B[K][J];
              }
            } else {
              if (NOUNIT) TEMP = TEMP * A[I][I].conjugate();
              for (K = I + 1; K <= M; K++) {
                TEMP = TEMP + A[K][I].conjugate() * B[K][J];
              }
            }
            B[I][J] = ALPHA * TEMP;
          }
        }
      }
    }
  } else {
    if (lsame(TRANSA, 'N')) {
      // Form  B := alpha*B*A.

      if (UPPER) {
        for (J = N; J >= 1; J--) {
          TEMP = ALPHA;
          if (NOUNIT) TEMP = TEMP * A[J][J];
          for (I = 1; I <= M; I++) {
            B[I][J] = TEMP * B[I][J];
          }
          for (K = 1; K <= J - 1; K++) {
            if (A[K][J] != Complex.zero) {
              TEMP = ALPHA * A[K][J];
              for (I = 1; I <= M; I++) {
                B[I][J] = B[I][J] + TEMP * B[I][K];
              }
            }
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          TEMP = ALPHA;
          if (NOUNIT) TEMP = TEMP * A[J][J];
          for (I = 1; I <= M; I++) {
            B[I][J] = TEMP * B[I][J];
          }
          for (K = J + 1; K <= N; K++) {
            if (A[K][J] != Complex.zero) {
              TEMP = ALPHA * A[K][J];
              for (I = 1; I <= M; I++) {
                B[I][J] = B[I][J] + TEMP * B[I][K];
              }
            }
          }
        }
      }
    } else {
      // Form  B := alpha*B*A**T   or   B := alpha*B*A**H.

      if (UPPER) {
        for (K = 1; K <= N; K++) {
          for (J = 1; J <= K - 1; J++) {
            if (A[J][K] != Complex.zero) {
              if (NOCONJ) {
                TEMP = ALPHA * A[J][K];
              } else {
                TEMP = ALPHA * A[J][K].conjugate();
              }
              for (I = 1; I <= M; I++) {
                B[I][J] = B[I][J] + TEMP * B[I][K];
              }
            }
          }
          TEMP = ALPHA;
          if (NOUNIT) {
            if (NOCONJ) {
              TEMP = TEMP * A[K][K];
            } else {
              TEMP = TEMP * A[K][K].conjugate();
            }
          }
          if (TEMP != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][K] = TEMP * B[I][K];
            }
          }
        }
      } else {
        for (K = N; K >= 1; K--) {
          for (J = K + 1; J <= N; J++) {
            if (A[J][K] != Complex.zero) {
              if (NOCONJ) {
                TEMP = ALPHA * A[J][K];
              } else {
                TEMP = ALPHA * A[J][K].conjugate();
              }
              for (I = 1; I <= M; I++) {
                B[I][J] = B[I][J] + TEMP * B[I][K];
              }
            }
          }
          TEMP = ALPHA;
          if (NOUNIT) {
            if (NOCONJ) {
              TEMP = TEMP * A[K][K];
            } else {
              TEMP = TEMP * A[K][K].conjugate();
            }
          }
          if (TEMP != Complex.one) {
            for (I = 1; I <= M; I++) {
              B[I][K] = TEMP * B[I][K];
            }
          }
        }
      }
    }
  }
}
