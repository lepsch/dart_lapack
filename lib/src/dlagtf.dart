// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlagtf(
  final int N,
  final Array<double> A_,
  final double LAMBDA,
  final Array<double> B_,
  final Array<double> C_,
  final double TOL,
  final Array<double> D_,
  final Array<int> IN_,
  final Box<int> INFO,
) {
  final A = A_.having();
  final B = B_.having();
  final C = C_.having();
  final D = D_.having();
  final IN = IN_.having();
  const ZERO = 0.0;
  int K;
  double EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL;

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('DLAGTF', -INFO.value);
    return;
  }

  if (N == 0) return;

  A[1] -= LAMBDA;
  IN[N] = 0;
  if (N == 1) {
    if (A[1] == ZERO) IN[1] = 1;
    return;
  }

  EPS = dlamch('Epsilon');

  TL = max(TOL, EPS);
  SCALE1 = A[1].abs() + B[1].abs();
  for (K = 1; K <= N - 1; K++) {
    A[K + 1] -= LAMBDA;
    SCALE2 = C[K].abs() + A[K + 1].abs();
    if (K < (N - 1)) SCALE2 += B[K + 1].abs();
    if (A[K] == ZERO) {
      PIV1 = ZERO;
    } else {
      PIV1 = A[K].abs() / SCALE1;
    }
    if (C[K] == ZERO) {
      IN[K] = 0;
      PIV2 = ZERO;
      SCALE1 = SCALE2;
      if (K < (N - 1)) D[K] = ZERO;
    } else {
      PIV2 = C[K].abs() / SCALE2;
      if (PIV2 <= PIV1) {
        IN[K] = 0;
        SCALE1 = SCALE2;
        C[K] /= A[K];
        A[K + 1] -= C[K] * B[K];
        if (K < (N - 1)) D[K] = ZERO;
      } else {
        IN[K] = 1;
        MULT = A[K] / C[K];
        A[K] = C[K];
        TEMP = A[K + 1];
        A[K + 1] = B[K] - MULT * TEMP;
        if (K < (N - 1)) {
          D[K] = B[K + 1];
          B[K + 1] = -MULT * D[K];
        }
        B[K] = TEMP;
        C[K] = MULT;
      }
    }
    if ((max(PIV1, PIV2) <= TL) && (IN[N] == 0)) IN[N] = K;
  }
  if ((A[N].abs() <= SCALE1 * TL) && (IN[N] == 0)) IN[N] = N;
}
