import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void dlatm5(
  final int PRTYPE,
  final int M,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final int LDB,
  final Matrix<double> C,
  final int LDC,
  final Matrix<double> D,
  final int LDD,
  final Matrix<double> E,
  final int LDE,
  final Matrix<double> F,
  final int LDF,
  final Matrix<double> R,
  final int LDR,
  final Matrix<double> L,
  final int LDL,
  final double ALPHA,
  final Box<int> QBLCKA,
  final Box<int> QBLCKB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0, TWENTY = 2.0e+1, HALF = 0.5, TWO = 2.0;
  int I, J, K;
  double IMEPS, REEPS;

  if (PRTYPE == 1) {
    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        if (I == J) {
          A[I][J] = ONE;
          D[I][J] = ONE;
        } else if (I == J - 1) {
          A[I][J] = -ONE;
          D[I][J] = ZERO;
        } else {
          A[I][J] = ZERO;
          D[I][J] = ZERO;
        }
      }
    }

    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        if (I == J) {
          B[I][J] = ONE - ALPHA;
          E[I][J] = ONE;
        } else if (I == J - 1) {
          B[I][J] = ONE;
          E[I][J] = ZERO;
        } else {
          B[I][J] = ZERO;
          E[I][J] = ZERO;
        }
      }
    }

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= N; J++) {
        R[I][J] = (HALF - sin((I / J).toDouble())) * TWENTY;
        L[I][J] = R[I][J];
      }
    }
  } else if (PRTYPE == 2 || PRTYPE == 3) {
    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        if (I <= J) {
          A[I][J] = (HALF - sin(I.toDouble())) * TWO;
          D[I][J] = (HALF - sin((I * J).toDouble())) * TWO;
        } else {
          A[I][J] = ZERO;
          D[I][J] = ZERO;
        }
      }
    }

    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        if (I <= J) {
          B[I][J] = (HALF - sin((I + J).toDouble())) * TWO;
          E[I][J] = (HALF - sin(J.toDouble())) * TWO;
        } else {
          B[I][J] = ZERO;
          E[I][J] = ZERO;
        }
      }
    }

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= N; J++) {
        R[I][J] = (HALF - sin((I * J).toDouble())) * TWENTY;
        L[I][J] = (HALF - sin((I + J).toDouble())) * TWENTY;
      }
    }

    if (PRTYPE == 3) {
      if (QBLCKA.value <= 1) QBLCKA.value = 2;
      for (K = 1;
          QBLCKA.value < 0 ? K >= M - 1 : K <= M - 1;
          K += QBLCKA.value) {
        A[K + 1][K + 1] = A[K][K];
        A[K + 1][K] = -sin(A[K][K + 1]);
      }

      if (QBLCKB.value <= 1) QBLCKB.value = 2;
      for (K = 1;
          QBLCKB.value < 0 ? K >= N - 1 : K <= N - 1;
          K += QBLCKB.value) {
        B[K + 1][K + 1] = B[K][K];
        B[K + 1][K] = -sin(B[K][K + 1]);
      }
    }
  } else if (PRTYPE == 4) {
    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        A[I][J] = (HALF - sin((I * J).toDouble())) * TWENTY;
        D[I][J] = (HALF - sin((I + J).toDouble())) * TWO;
      }
    }

    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        B[I][J] = (HALF - sin((I + J).toDouble())) * TWENTY;
        E[I][J] = (HALF - sin((I * J).toDouble())) * TWO;
      }
    }

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= N; J++) {
        R[I][J] = (HALF - sin((J / I).toDouble())) * TWENTY;
        L[I][J] = (HALF - sin((I * J).toDouble())) * TWO;
      }
    }
  } else if (PRTYPE >= 5) {
    REEPS = HALF * TWO * TWENTY / ALPHA;
    IMEPS = (HALF - TWO) / ALPHA;
    for (I = 1; I <= M; I++) {
      for (J = 1; J <= N; J++) {
        R[I][J] = (HALF - sin((I * J).toDouble())) * ALPHA / TWENTY;
        L[I][J] = (HALF - sin((I + J).toDouble())) * ALPHA / TWENTY;
      }
    }

    for (I = 1; I <= M; I++) {
      D[I][I] = ONE;
    }

    for (I = 1; I <= M; I++) {
      if (I <= 4) {
        A[I][I] = ONE;
        if (I > 2) A[I][I] = ONE + REEPS;
        if ((I % 2) != 0 && I < M) {
          A[I][I + 1] = IMEPS;
        } else if (I > 1) {
          A[I][I - 1] = -IMEPS;
        }
      } else if (I <= 8) {
        if (I <= 6) {
          A[I][I] = REEPS;
        } else {
          A[I][I] = -REEPS;
        }
        if ((I % 2) != 0 && I < M) {
          A[I][I + 1] = ONE;
        } else if (I > 1) {
          A[I][I - 1] = -ONE;
        }
      } else {
        A[I][I] = ONE;
        if ((I % 2) != 0 && I < M) {
          A[I][I + 1] = IMEPS * 2;
        } else if (I > 1) {
          A[I][I - 1] = -IMEPS * 2;
        }
      }
    }

    for (I = 1; I <= N; I++) {
      E[I][I] = ONE;
      if (I <= 4) {
        B[I][I] = -ONE;
        if (I > 2) B[I][I] = ONE - REEPS;
        if ((I % 2) != 0 && I < N) {
          B[I][I + 1] = IMEPS;
        } else if (I > 1) {
          B[I][I - 1] = -IMEPS;
        }
      } else if (I <= 8) {
        if (I <= 6) {
          B[I][I] = REEPS;
        } else {
          B[I][I] = -REEPS;
        }
        if ((I % 2) != 0 && I < N) {
          B[I][I + 1] = ONE + IMEPS;
        } else if (I > 1) {
          B[I][I - 1] = -ONE - IMEPS;
        }
      } else {
        B[I][I] = ONE - REEPS;
        if ((I % 2) != 0 && I < N) {
          B[I][I + 1] = IMEPS * 2;
        } else if (I > 1) {
          B[I][I - 1] = -IMEPS * 2;
        }
      }
    }
  }

  // Compute rhs (C, F)

  dgemm('N', 'N', M, N, M, ONE, A, LDA, R, LDR, ZERO, C, LDC);
  dgemm('N', 'N', M, N, N, -ONE, L, LDL, B, LDB, ONE, C, LDC);
  dgemm('N', 'N', M, N, M, ONE, D, LDD, R, LDR, ZERO, F, LDF);
  dgemm('N', 'N', M, N, N, -ONE, L, LDL, E, LDE, ONE, F, LDF);
}
