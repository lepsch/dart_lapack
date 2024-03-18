import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void zlaptm(
  final String UPLO,
  final int N,
  final int NRHS,
  final double ALPHA,
  final Array<double> D_,
  final Array<Complex> E_,
  final Matrix<Complex> X_,
  final int LDX,
  final double BETA,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, ZERO = 0.0;

  if (N == 0) return;

  if (BETA == ZERO) {
    for (var J = 1; J <= NRHS; J++) {
      for (var I = 1; I <= N; I++) {
        B[I][J] = Complex.zero;
      }
    }
  } else if (BETA == -ONE) {
    for (var J = 1; J <= NRHS; J++) {
      for (var I = 1; I <= N; I++) {
        B[I][J] = -B[I][J];
      }
    }
  }

  if (ALPHA == ONE) {
    if (lsame(UPLO, 'U')) {
      // Compute B := B + A*X, where E is the superdiagonal of A.

      for (var J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] += D[1].toComplex() * X[1][J];
        } else {
          B[1][J] += D[1].toComplex() * X[1][J] + E[1] * X[2][J];
          B[N][J] +=
              E[N - 1].conjugate() * X[N - 1][J] + D[N].toComplex() * X[N][J];
          for (var I = 2; I <= N - 1; I++) {
            B[I][J] += E[I - 1].conjugate() * X[I - 1][J] +
                D[I].toComplex() * X[I][J] +
                E[I] * X[I + 1][J];
          }
        }
      }
    } else {
      // Compute B := B + A*X, where E is the subdiagonal of A.

      for (var J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] += D[1].toComplex() * X[1][J];
        } else {
          B[1][J] += D[1].toComplex() * X[1][J] + E[1].conjugate() * X[2][J];
          B[N][J] += E[N - 1] * X[N - 1][J] + D[N].toComplex() * X[N][J];
          for (var I = 2; I <= N - 1; I++) {
            B[I][J] += E[I - 1] * X[I - 1][J] +
                D[I].toComplex() * X[I][J] +
                E[I].conjugate() * X[I + 1][J];
          }
        }
      }
    }
  } else if (ALPHA == -ONE) {
    if (lsame(UPLO, 'U')) {
      // Compute B := B - A*X, where E is the superdiagonal of A.

      for (var J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] -= D[1].toComplex() * X[1][J];
        } else {
          B[1][J] -= D[1].toComplex() * X[1][J] + E[1] * X[2][J];
          B[N][J] -=
              E[N - 1].conjugate() * X[N - 1][J] + D[N].toComplex() * X[N][J];
          for (var I = 2; I <= N - 1; I++) {
            B[I][J] -= E[I - 1].conjugate() * X[I - 1][J] +
                D[I].toComplex() * X[I][J] +
                E[I] * X[I + 1][J];
          }
        }
      }
    } else {
      // Compute B := B - A*X, where E is the subdiagonal of A.

      for (var J = 1; J <= NRHS; J++) {
        if (N == 1) {
          B[1][J] -= D[1].toComplex() * X[1][J];
        } else {
          B[1][J] -= D[1].toComplex() * X[1][J] + E[1].conjugate() * X[2][J];
          B[N][J] -= E[N - 1] * X[N - 1][J] + D[N].toComplex() * X[N][J];
          for (var I = 2; I <= N - 1; I++) {
            B[I][J] -= E[I - 1] * X[I - 1][J] +
                D[I].toComplex() * X[I][J] +
                E[I].conjugate() * X[I + 1][J];
          }
        }
      }
    }
  }
}
