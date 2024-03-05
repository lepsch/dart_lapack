import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlatm5(
  final int PRTYPE,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> D_,
  final int LDD,
  final Matrix<Complex> E_,
  final int LDE,
  final Matrix<Complex> F_,
  final int LDF,
  final Matrix<Complex> R_,
  final int LDR,
  final Matrix<Complex> L_,
  final int LDL,
  final double ALPHA,
  final Box<int> QBLCKA,
  final Box<int> QBLCKB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final D = D_.having(ld: LDD);
  final E = E_.having(ld: LDE);
  final F = F_.having(ld: LDF);
  final R = R_.having(ld: LDR);
  final L = L_.having(ld: LDL);
  const TWO = Complex(2.0, 0.0),
      HALF = Complex(0.5, 0.0),
      TWENTY = Complex(2.0e+1, 0.0);
  int I, J, K;
  Complex IMEPS, REEPS;

  if (PRTYPE == 1) {
    for (I = 1; I <= M; I++) {
      // 20
      for (J = 1; J <= M; J++) {
        // 10
        if (I == J) {
          A[I][J] = Complex.one;
          D[I][J] = Complex.one;
        } else if (I == J - 1) {
          A[I][J] = -Complex.one;
          D[I][J] = Complex.zero;
        } else {
          A[I][J] = Complex.zero;
          D[I][J] = Complex.zero;
        }
      } // 10
    } // 20

    for (I = 1; I <= N; I++) {
      // 40
      for (J = 1; J <= N; J++) {
        // 30
        if (I == J) {
          B[I][J] = Complex.one - ALPHA.toComplex();
          E[I][J] = Complex.one;
        } else if (I == J - 1) {
          B[I][J] = Complex.one;
          E[I][J] = Complex.zero;
        } else {
          B[I][J] = Complex.zero;
          E[I][J] = Complex.zero;
        }
      } // 30
    } // 40

    for (I = 1; I <= M; I++) {
      // 60
      for (J = 1; J <= N; J++) {
        // 50
        R[I][J] = (HALF - (I / J).toComplex()).sin() * TWENTY;
        L[I][J] = R[I][J];
      } // 50
    } // 60
  } else if (PRTYPE == 2 || PRTYPE == 3) {
    for (I = 1; I <= M; I++) {
      // 80
      for (J = 1; J <= M; J++) {
        // 70
        if (I <= J) {
          A[I][J] = (HALF - I.toComplex().sin()) * TWO;
          D[I][J] = (HALF - (I * J).toComplex().sin()) * TWO;
        } else {
          A[I][J] = Complex.zero;
          D[I][J] = Complex.zero;
        }
      } // 70
    } // 80

    for (I = 1; I <= N; I++) {
      // 100
      for (J = 1; J <= N; J++) {
        // 90
        if (I <= J) {
          B[I][J] = (HALF - (I + J).toComplex()).sin() * TWO;
          E[I][J] = (HALF - J.toComplex()).sin() * TWO;
        } else {
          B[I][J] = Complex.zero;
          E[I][J] = Complex.zero;
        }
      } // 90
    } // 100

    for (I = 1; I <= M; I++) {
      // 120
      for (J = 1; J <= N; J++) {
        // 110
        R[I][J] = (HALF - (I * J).toComplex().sin()) * TWENTY;
        L[I][J] = (HALF - (I + J).toComplex().sin()) * TWENTY;
      } // 110
    } // 120

    if (PRTYPE == 3) {
      if (QBLCKA.value <= 1) QBLCKA.value = 2;
      for (K = 1;
          QBLCKA.value < 0 ? K >= M - 1 : K <= M - 1;
          K += QBLCKA.value) {
        // 130
        A[K + 1][K + 1] = A[K][K];
        A[K + 1][K] = -A[K][K + 1].sin();
      } // 130

      if (QBLCKB.value <= 1) QBLCKB.value = 2;
      for (K = 1;
          QBLCKB.value < 0 ? K >= N - 1 : K <= N - 1;
          K += QBLCKB.value) {
        // 140
        B[K + 1][K + 1] = B[K][K];
        B[K + 1][K] = -B[K][K + 1].sin();
      } // 140
    }
  } else if (PRTYPE == 4) {
    for (I = 1; I <= M; I++) {
      // 160
      for (J = 1; J <= M; J++) {
        // 150
        A[I][J] = (HALF - (I * J).toComplex().sin()) * TWENTY;
        D[I][J] = (HALF - (I + J).toComplex().sin()) * TWO;
      } // 150
    } // 160

    for (I = 1; I <= N; I++) {
      // 180
      for (J = 1; J <= N; J++) {
        // 170
        B[I][J] = (HALF - (I + J).toComplex().sin()) * TWENTY;
        E[I][J] = (HALF - (I * J).toComplex().sin()) * TWO;
      } // 170
    } // 180

    for (I = 1; I <= M; I++) {
      // 200
      for (J = 1; J <= N; J++) {
        // 190
        R[I][J] = (HALF - (J / I).toComplex().sin()) * TWENTY;
        L[I][J] = (HALF - (I * J).toComplex().sin()) * TWO;
      } // 190
    } // 200
  } else if (PRTYPE >= 5) {
    REEPS = HALF * TWO * TWENTY / ALPHA.toComplex();
    IMEPS = (HALF - TWO) / ALPHA.toComplex();
    for (I = 1; I <= M; I++) {
      // 220
      for (J = 1; J <= N; J++) {
        // 210
        R[I][J] =
            (HALF - (I * J).toComplex().sin()) * ALPHA.toComplex() / TWENTY;
        L[I][J] =
            (HALF - (I + J).toComplex().sin()) * ALPHA.toComplex() / TWENTY;
      } // 210
    } // 220

    for (I = 1; I <= M; I++) {
      // 230
      D[I][I] = Complex.one;
    } // 230

    for (I = 1; I <= M; I++) {
      // 240
      if (I <= 4) {
        A[I][I] = Complex.one;
        if (I > 2) A[I][I] = Complex.one + REEPS;
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
          A[I][I + 1] = Complex.one;
        } else if (I > 1) {
          A[I][I - 1] = -Complex.one;
        }
      } else {
        A[I][I] = Complex.one;
        if ((I % 2) != 0 && I < M) {
          A[I][I + 1] = IMEPS * 2.toComplex();
        } else if (I > 1) {
          A[I][I - 1] = -IMEPS * 2.toComplex();
        }
      }
    } // 240

    for (I = 1; I <= N; I++) {
      // 250
      E[I][I] = Complex.one;
      if (I <= 4) {
        B[I][I] = -Complex.one;
        if (I > 2) B[I][I] = Complex.one - REEPS;
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
          B[I][I + 1] = Complex.one + IMEPS;
        } else if (I > 1) {
          B[I][I - 1] = -Complex.one - IMEPS;
        }
      } else {
        B[I][I] = Complex.one - REEPS;
        if ((I % 2) != 0 && I < N) {
          B[I][I + 1] = IMEPS * 2.toComplex();
        } else if (I > 1) {
          B[I][I - 1] = -IMEPS * 2.toComplex();
        }
      }
    } // 250
  }

  // Compute rhs (C, F)

  zgemm('N', 'N', M, N, M, Complex.one, A, LDA, R, LDR, Complex.zero, C, LDC);
  zgemm('N', 'N', M, N, N, -Complex.one, L, LDL, B, LDB, Complex.one, C, LDC);
  zgemm('N', 'N', M, N, M, Complex.one, D, LDD, R, LDR, Complex.zero, F, LDF);
  zgemm('N', 'N', M, N, N, -Complex.one, L, LDL, E, LDE, Complex.one, F, LDF);
}
