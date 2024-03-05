import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

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
  double TEMP1, TEMP2;
  int I, INFO, J, K, NROWA;
  bool UPPER;
  const ONE = 1.0, ZERO = 0.0;

  // Set NROWA as the number of rows of A.

  if (lsame(SIDE, 'L')) {
    NROWA = M;
  } else {
    NROWA = N;
  }
  UPPER = lsame(UPLO, 'U');

  // Test the input parameters.

  INFO = 0;
  if ((!lsame(SIDE, 'L')) && (!lsame(SIDE, 'R'))) {
    INFO = 1;
  } else if ((!UPPER) && (!lsame(UPLO, 'L'))) {
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
    xerbla('DSYMM ', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // And when  alpha == zero.

  if (ALPHA == ZERO) {
    if (BETA == ZERO) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          C[I][J] = ZERO;
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

  if (lsame(SIDE, 'L')) {
    // Form  C := alpha*A*B + beta*C.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP1 = ALPHA * B[I][J];
          TEMP2 = ZERO;
          for (K = 1; K <= I - 1; K++) {
            C[K][J] = C[K][J] + TEMP1 * A[K][I];
            TEMP2 = TEMP2 + B[K][J] * A[K][I];
          }
          if (BETA == ZERO) {
            C[I][J] = TEMP1 * A[I][I] + ALPHA * TEMP2;
          } else {
            C[I][J] = BETA * C[I][J] + TEMP1 * A[I][I] + ALPHA * TEMP2;
          }
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        for (I = M; I >= 1; I--) {
          TEMP1 = ALPHA * B[I][J];
          TEMP2 = ZERO;
          for (K = I + 1; K <= M; K++) {
            C[K][J] = C[K][J] + TEMP1 * A[K][I];
            TEMP2 = TEMP2 + B[K][J] * A[K][I];
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

    for (J = 1; J <= N; J++) {
      TEMP1 = ALPHA * A[J][J];
      if (BETA == ZERO) {
        for (I = 1; I <= M; I++) {
          C[I][J] = TEMP1 * B[I][J];
        }
      } else {
        for (I = 1; I <= M; I++) {
          C[I][J] = BETA * C[I][J] + TEMP1 * B[I][J];
        }
      }
      for (K = 1; K <= J - 1; K++) {
        if (UPPER) {
          TEMP1 = ALPHA * A[K][J];
        } else {
          TEMP1 = ALPHA * A[J][K];
        }
        for (I = 1; I <= M; I++) {
          C[I][J] = C[I][J] + TEMP1 * B[I][K];
        }
      }
      for (K = J + 1; K <= N; K++) {
        if (UPPER) {
          TEMP1 = ALPHA * A[J][K];
        } else {
          TEMP1 = ALPHA * A[K][J];
        }
        for (I = 1; I <= M; I++) {
          C[I][J] = C[I][J] + TEMP1 * B[I][K];
        }
      }
    }
  }
}
