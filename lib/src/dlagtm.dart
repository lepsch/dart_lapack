import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlagtm(
  final String TRANS,
  final int N,
  final int NRHS,
  final double ALPHA,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Matrix<double> X_,
  final int LDX,
  final double BETA,
  final Matrix<double> B_,
  final int LDB,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.dim();
  final D = D_.dim();
  final DU = DU_.dim();
  final X = X_.dim(LDX);
  final B = B_.dim(LDB);
  const ONE = 1.0, ZERO = 0.0;
  int I, J;

  if (N == 0) return;

  // Multiply B by BETA if BETA != 1.

  if (BETA == ZERO) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        B[I][J] = ZERO;
      }
    }
  } else if (BETA == -ONE) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        B[I][J] = -B[I][J];
      }
    }
  }

  if (ALPHA == ONE) {
    if (lsame(TRANS, 'N')) {
      // Compute B := B + A*X

      for (J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] = B[1][J] + D[1] * X[1][J];
        } else {
          B[1][J] = B[1][J] + D[1] * X[1][J] + DU[1] * X[2][J];
          B[N][J] = B[N][J] + DL[N - 1] * X[N - 1][J] + D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            B[I][J] = B[I][J] +
                DL[I - 1] * X[I - 1][J] +
                D[I] * X[I][J] +
                DU[I] * X[I + 1][J];
          }
        }
      }
    } else {
      // Compute B := B + A**T*X

      for (J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] = B[1][J] + D[1] * X[1][J];
        } else {
          B[1][J] = B[1][J] + D[1] * X[1][J] + DL[1] * X[2][J];
          B[N][J] = B[N][J] + DU[N - 1] * X[N - 1][J] + D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            B[I][J] = B[I][J] +
                DU[I - 1] * X[I - 1][J] +
                D[I] * X[I][J] +
                DL[I] * X[I + 1][J];
          }
        }
      }
    }
  } else if (ALPHA == -ONE) {
    if (lsame(TRANS, 'N')) {
      // Compute B := B - A*X

      for (J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] = B[1][J] - D[1] * X[1][J];
        } else {
          B[1][J] = B[1][J] - D[1] * X[1][J] - DU[1] * X[2][J];
          B[N][J] = B[N][J] - DL[N - 1] * X[N - 1][J] - D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            B[I][J] = B[I][J] -
                DL[I - 1] * X[I - 1][J] -
                D[I] * X[I][J] -
                DU[I] * X[I + 1][J];
          }
        }
      }
    } else {
      // Compute B := B - A**T*X

      for (J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] = B[1][J] - D[1] * X[1][J];
        } else {
          B[1][J] = B[1][J] - D[1] * X[1][J] - DL[1] * X[2][J];
          B[N][J] = B[N][J] - DU[N - 1] * X[N - 1][J] - D[N] * X[N][J];
          for (I = 2; I <= N - 1; I++) {
            B[I][J] = B[I][J] -
                DU[I - 1] * X[I - 1][J] -
                D[I] * X[I][J] -
                DL[I] * X[I + 1][J];
          }
        }
      }
    }
  }
}
