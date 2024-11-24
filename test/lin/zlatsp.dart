// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

import '../matgen/zlarnd.dart';

void zlatsp(
  final String UPLO,
  final int N,
  final Array<Complex> X_,
  final Array<int> ISEED_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final ISEED = ISEED_.having();
  const EYE = Complex(0.0, 1.0), TWO = Complex(2.0);

  // Initialize constants

  final ALPHA = ((1.0 + sqrt(17.0)) / 8.0).toComplex();
  final BETA = ALPHA - (1.0 / 1000.0).toComplex();
  final ALPHA3 = ALPHA * ALPHA * ALPHA;

  // Fill the matrix with zeros.

  for (var J = 1; J <= N * (N + 1) ~/ 2; J++) {
    X[J] = Complex.zero;
  }

  // UPLO = 'U':  Upper triangular storage

  if (UPLO == 'U') {
    final N5 = N - 5 * (N ~/ 5) + 1;

    var JJ = N * (N + 1) ~/ 2;
    for (var J = N; J >= N5; J -= 5) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[JJ] = A;
      X[JJ - 2] = B;
      JJ -= J;
      X[JJ] = zlarnd(2, ISEED);
      X[JJ - 1] = R;
      JJ -= (J - 1);
      X[JJ] = C;
      JJ -= (J - 2);
      X[JJ] = zlarnd(2, ISEED);
      JJ -= (J - 3);
      X[JJ] = zlarnd(2, ISEED);
      if (X[JJ + (J - 3)].abs() > X[JJ].abs()) {
        X[JJ + (J - 4)] = TWO * X[JJ + (J - 3)];
      } else {
        X[JJ + (J - 4)] = TWO * X[JJ];
      }
      JJ -= (J - 4);
    }

    // Clean-up for N not a multiple of 5.

    var J = N5 - 1;
    if (J > 2) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[JJ] = A;
      X[JJ - 2] = B;
      JJ -= J;
      X[JJ] = zlarnd(2, ISEED);
      X[JJ - 1] = R;
      JJ -= (J - 1);
      X[JJ] = C;
      JJ -= (J - 2);
      J -= 3;
    }
    if (J > 1) {
      X[JJ] = zlarnd(2, ISEED);
      X[JJ - J] = zlarnd(2, ISEED);
      if (X[JJ].abs() > X[JJ - J].abs()) {
        X[JJ - 1] = TWO * X[JJ];
      } else {
        X[JJ - 1] = TWO * X[JJ - J];
      }
      JJ -= J + (J - 1);
      J -= 2;
    } else if (J == 1) {
      X[JJ] = zlarnd(2, ISEED);
      J--;
    }

    // UPLO = 'L':  Lower triangular storage
  } else {
    final N5 = (N ~/ 5) * 5;

    var JJ = 1;
    for (var J = 1; J <= N5; J += 5) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[JJ] = A;
      X[JJ + 2] = B;
      JJ += (N - J + 1);
      X[JJ] = zlarnd(2, ISEED);
      X[JJ + 1] = R;
      JJ += (N - J);
      X[JJ] = C;
      JJ += (N - J - 1);
      X[JJ] = zlarnd(2, ISEED);
      JJ += (N - J - 2);
      X[JJ] = zlarnd(2, ISEED);
      if (X[JJ - (N - J - 2)].abs() > X[JJ].abs()) {
        X[JJ - (N - J - 2) + 1] = TWO * X[JJ - (N - J - 2)];
      } else {
        X[JJ - (N - J - 2) + 1] = TWO * X[JJ];
      }
      JJ += (N - J - 3);
    }

    // Clean-up for N not a multiple of 5.

    var J = N5 + 1;
    if (J < N - 1) {
      final A = ALPHA3 * zlarnd(5, ISEED);
      final B = zlarnd(5, ISEED) / ALPHA;
      final C = A - TWO * B * EYE;
      final R = C / BETA;
      X[JJ] = A;
      X[JJ + 2] = B;
      JJ += (N - J + 1);
      X[JJ] = zlarnd(2, ISEED);
      X[JJ + 1] = R;
      JJ += (N - J);
      X[JJ] = C;
      JJ += (N - J - 1);
      J += 3;
    }
    if (J < N) {
      X[JJ] = zlarnd(2, ISEED);
      X[JJ + (N - J + 1)] = zlarnd(2, ISEED);
      if (X[JJ].abs() > X[JJ + (N - J + 1)].abs()) {
        X[JJ + 1] = TWO * X[JJ];
      } else {
        X[JJ + 1] = TWO * X[JJ + (N - J + 1)];
      }
      JJ += (N - J + 1) + (N - J);
      J += 2;
    } else if (J == N) {
      X[JJ] = zlarnd(2, ISEED);
      JJ += (N - J + 1);
      J++;
    }
  }
}
