// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlaqz1(
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final double SR1,
  final double SR2,
  final double SI,
  final double BETA1,
  final double BETA2,
  final Array<double> V,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ZERO = 0.0, ONE = 1.0;
  double SAFMIN, SAFMAX, SCALE1, SCALE2;
  final W = Array<double>(2);

  SAFMIN = dlamch('SAFE MINIMUM');
  SAFMAX = ONE / SAFMIN;

  // Calculate first shifted vector

  W[1] = BETA1 * A[1][1] - SR1 * B[1][1];
  W[2] = BETA1 * A[2][1] - SR1 * B[2][1];
  SCALE1 = sqrt(W[1].abs()) * sqrt(W[2].abs());
  if (SCALE1 >= SAFMIN && SCALE1 <= SAFMAX) {
    W[1] /= SCALE1;
    W[2] /= SCALE1;
  }

  // Solve linear system

  W[2] /= B[2][2];
  W[1] = (W[1] - B[1][2] * W[2]) / B[1][1];
  SCALE2 = sqrt(W[1].abs()) * sqrt(W[2].abs());
  if (SCALE2 >= SAFMIN && SCALE2 <= SAFMAX) {
    W[1] /= SCALE2;
    W[2] /= SCALE2;
  }

  // Apply second shift

  V[1] = BETA2 * (A[1][1] * W[1] + A[1][2] * W[2]) -
      SR2 * (B[1][1] * W[1] + B[1][2] * W[2]);
  V[2] = BETA2 * (A[2][1] * W[1] + A[2][2] * W[2]) -
      SR2 * (B[2][1] * W[1] + B[2][2] * W[2]);

  V[3] = BETA2 * (A[3][1] * W[1] + A[3][2] * W[2]) -
      SR2 * (B[3][1] * W[1] + B[3][2] * W[2]);

  // Account for imaginary part

  V[1] += SI * SI * B[1][1] / SCALE1 / SCALE2;

  // Check for overflow

  if (V[1].abs() > SAFMAX ||
      V[2].abs() > SAFMAX ||
      V[3].abs() > SAFMAX ||
      disnan(V[1]) ||
      disnan(V[2]) ||
      disnan(V[3])) {
    V[1] = ZERO;
    V[2] = ZERO;
    V[3] = ZERO;
  }
}
