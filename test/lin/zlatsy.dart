import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

import '../matgen/zlarnd.dart';

void zlatsy(
  final String UPLO,
  final int N,
  final Matrix<Complex> X_,
  final int LDX,
  final Array<int> ISEED_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having(ld: LDX);
  final ISEED = ISEED_.having();
  const EYE = Complex(0.0, 1.0), TWO = Complex(2.0);

  // Initialize constants

  final ALPHA = ((1.0 + sqrt(17.0)) / 8.0).toComplex();
  final BETA = ALPHA - (1.0 / 1000.0).toComplex();
  final ALPHA3 = ALPHA * ALPHA * ALPHA;

  // UPLO = 'U':  Upper triangular storage

  if (UPLO == 'U') {
    // Fill the upper triangle of the matrix with zeros.

    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J; I++) {
        X[I][J] = Complex.zero;
      }
    }
    final N5 = N - 5 * (N ~/ 5) + 1;

    for (var I = N; I >= N5; I -= 5) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[I][I] = A;
      X[I - 2][I] = B;
      X[I - 2][I - 1] = R;
      X[I - 2][I - 2] = C;
      X[I - 1][I - 1] = zlarnd(2, ISEED);
      X[I - 3][I - 3] = zlarnd(2, ISEED);
      X[I - 4][I - 4] = zlarnd(2, ISEED);
      if (X[I - 3][I - 3].abs() > X[I - 4][I - 4].abs()) {
        X[I - 4][I - 3] = TWO * X[I - 3][I - 3];
      } else {
        X[I - 4][I - 3] = TWO * X[I - 4][I - 4];
      }
    }

    // Clean-up for N not a multiple of 5.

    var I = N5 - 1;
    if (I > 2) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[I][I] = A;
      X[I - 2][I] = B;
      X[I - 2][I - 1] = R;
      X[I - 2][I - 2] = C;
      X[I - 1][I - 1] = zlarnd(2, ISEED);
      I = I - 3;
    }
    if (I > 1) {
      X[I][I] = zlarnd(2, ISEED);
      X[I - 1][I - 1] = zlarnd(2, ISEED);
      if (X[I][I].abs() > X[I - 1][I - 1].abs()) {
        X[I - 1][I] = TWO * X[I][I];
      } else {
        X[I - 1][I] = TWO * X[I - 1][I - 1];
      }
      I = I - 2;
    } else if (I == 1) {
      X[I][I] = zlarnd(2, ISEED);
      I--;
    }

    // UPLO = 'L':  Lower triangular storage
  } else {
    // Fill the lower triangle of the matrix with zeros.

    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= N; I++) {
        X[I][J] = Complex.zero;
      }
    }
    final N5 = (N ~/ 5) * 5;

    for (var I = 1; I <= N5; I += 5) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[I][I] = A;
      X[I + 2][I] = B;
      X[I + 2][I + 1] = R;
      X[I + 2][I + 2] = C;
      X[I + 1][I + 1] = zlarnd(2, ISEED);
      X[I + 3][I + 3] = zlarnd(2, ISEED);
      X[I + 4][I + 4] = zlarnd(2, ISEED);
      if (X[I + 3][I + 3].abs() > X[I + 4][I + 4].abs()) {
        X[I + 4][I + 3] = TWO * X[I + 3][I + 3];
      } else {
        X[I + 4][I + 3] = TWO * X[I + 4][I + 4];
      }
    }

    // Clean-up for N not a multiple of 5.

    var I = N5 + 1;
    if (I < N - 1) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[I][I] = A;
      X[I + 2][I] = B;
      X[I + 2][I + 1] = R;
      X[I + 2][I + 2] = C;
      X[I + 1][I + 1] = zlarnd(2, ISEED);
      I = I + 3;
    }
    if (I < N) {
      X[I][I] = zlarnd(2, ISEED);
      X[I + 1][I + 1] = zlarnd(2, ISEED);
      if (X[I][I].abs() > X[I + 1][I + 1].abs()) {
        X[I + 1][I] = TWO * X[I][I];
      } else {
        X[I + 1][I] = TWO * X[I + 1][I + 1];
      }
      I = I + 2;
    } else if (I == N) {
      X[I][I] = zlarnd(2, ISEED);
      I++;
    }
  }
}
