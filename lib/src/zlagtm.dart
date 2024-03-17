import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlagtm(
  final String TRANS,
  final int N,
  final int NRHS,
  final double ALPHA,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Matrix<Complex> X_,
  final int LDX,
  final double BETA,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, ZERO = 0.0;
  int I, J;

  if (N == 0) return;

  // Multiply B by BETA if BETA != 1.

  if (BETA == ZERO) {
    for (J = 1; J <= NRHS; J++) {
      // 20
      for (I = 1; I <= N; I++) {
        // 10
        B[I][J] = Complex.zero;
      } // 10
    } // 20
  } else if (BETA == -ONE) {
    for (J = 1; J <= NRHS; J++) {
      // 40
      for (I = 1; I <= N; I++) {
        // 30
        B[I][J] = -B[I][J];
      } // 30
    } // 40
  }

  if (ALPHA == ONE) {
    if (lsame(TRANS, 'N')) {
      // Compute B := B + A*X

      for (J = 1; J <= NRHS; J++) {
        // 60
        if (N == 1) {
          B[1][J] += D[1] * X[1][J];
        } else {
          B[1][J] += D[1] * X[1][J] + DU[1] * X[2][J];
          B[N][J] += DL[N - 1] * X[N - 1][J] + D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            // 50
            B[I][J] +=
                DL[I - 1] * X[I - 1][J] + D[I] * X[I][J] + DU[I] * X[I + 1][J];
          } // 50
        }
      } // 60
    } else if (lsame(TRANS, 'T')) {
      // Compute B := B + A**T * X

      for (J = 1; J <= NRHS; J++) {
        // 80
        if (N == 1) {
          B[1][J] += D[1] * X[1][J];
        } else {
          B[1][J] += D[1] * X[1][J] + DL[1] * X[2][J];
          B[N][J] += DU[N - 1] * X[N - 1][J] + D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            // 70
            B[I][J] +=
                DU[I - 1] * X[I - 1][J] + D[I] * X[I][J] + DL[I] * X[I + 1][J];
          } // 70
        }
      } // 80
    } else if (lsame(TRANS, 'C')) {
      // Compute B := B + A**H * X

      for (J = 1; J <= NRHS; J++) {
        // 100
        if (N == 1) {
          B[1][J] += (D[1].conjugate()) * X[1][J];
        } else {
          B[1][J] +=
              (D[1].conjugate()) * X[1][J] + (DL[1].conjugate()) * X[2][J];
          B[N][J] += (DU[N - 1].conjugate()) * X[N - 1][J] +
              (D[N].conjugate()) * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            // 90
            B[I][J] += (DU[I - 1].conjugate()) * X[I - 1][J] +
                (D[I].conjugate()) * X[I][J] +
                (DL[I].conjugate()) * X[I + 1][J];
          } // 90
        }
      } // 100
    }
  } else if (ALPHA == -ONE) {
    if (lsame(TRANS, 'N')) {
      // Compute B := B - A*X

      for (J = 1; J <= NRHS; J++) {
        // 120
        if (N == 1) {
          B[1][J] -= D[1] * X[1][J];
        } else {
          B[1][J] -= D[1] * X[1][J] - DU[1] * X[2][J];
          B[N][J] -= DL[N - 1] * X[N - 1][J] - D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            // 110
            B[I][J] -=
                DL[I - 1] * X[I - 1][J] - D[I] * X[I][J] - DU[I] * X[I + 1][J];
          } // 110
        }
      } // 120
    } else if (lsame(TRANS, 'T')) {
      // Compute B := B - A**T *X

      for (J = 1; J <= NRHS; J++) {
        // 140
        if (N == 1) {
          B[1][J] -= D[1] * X[1][J];
        } else {
          B[1][J] -= D[1] * X[1][J] - DL[1] * X[2][J];
          B[N][J] -= DU[N - 1] * X[N - 1][J] - D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            // 130
            B[I][J] -=
                DU[I - 1] * X[I - 1][J] - D[I] * X[I][J] - DL[I] * X[I + 1][J];
          } // 130
        }
      } // 140
    } else if (lsame(TRANS, 'C')) {
      // Compute B := B - A**H *X

      for (J = 1; J <= NRHS; J++) {
        // 160
        if (N == 1) {
          B[1][J] -= (D[1].conjugate()) * X[1][J];
        } else {
          B[1][J] -=
              (D[1].conjugate()) * X[1][J] - (DL[1].conjugate()) * X[2][J];
          B[N][J] -= (DU[N - 1].conjugate()) * X[N - 1][J] -
              (D[N].conjugate()) * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            // 150
            B[I][J] -= (DU[I - 1].conjugate()) * X[I - 1][J] -
                (D[I].conjugate()) * X[I][J] -
                (DL[I].conjugate()) * X[I + 1][J];
          } // 150
        }
      } // 160
    }
  }
}
