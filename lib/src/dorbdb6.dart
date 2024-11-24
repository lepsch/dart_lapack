// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlassq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorbdb6(
  final int M1,
  final int M2,
  final int N,
  final Array<double> X1_,
  final int INCX1,
  final Array<double> X2_,
  final int INCX2,
  final Matrix<double> Q1_,
  final int LDQ1,
  final Matrix<double> Q2_,
  final int LDQ2,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X1 = X1_.having();
  final X2 = X2_.having();
  final Q1 = Q1_.having(ld: LDQ1);
  final Q2 = Q2_.having(ld: LDQ2);
  final WORK = WORK_.having();
  const ALPHA = 0.83,
      // REALONE = 1.0,
      REALZERO = 0.0;
  const NEGONE = -1.0, ONE = 1.0, ZERO = 0.0;
  int I, IX;
  double EPS, NORM, NORM_NEW;
  final SCL = Box(0.0), SSQ = Box(0.0);

  // Test input arguments

  INFO.value = 0;
  if (M1 < 0) {
    INFO.value = -1;
  } else if (M2 < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (INCX1 < 1) {
    INFO.value = -5;
  } else if (INCX2 < 1) {
    INFO.value = -7;
  } else if (LDQ1 < max(1, M1)) {
    INFO.value = -9;
  } else if (LDQ2 < max(1, M2)) {
    INFO.value = -11;
  } else if (LWORK < N) {
    INFO.value = -13;
  }

  if (INFO.value != 0) {
    xerbla('DORBDB6', -INFO.value);
    return;
  }

  EPS = dlamch('Precision');

  // Compute the Euclidean norm of X

  SCL.value = REALZERO;
  SSQ.value = REALZERO;
  dlassq(M1, X1, INCX1, SCL, SSQ);
  dlassq(M2, X2, INCX2, SCL, SSQ);
  NORM = SCL.value * sqrt(SSQ.value);

  // First, project X onto the orthogonal complement of Q's column
  // space

  if (M1 == 0) {
    for (I = 1; I <= N; I++) {
      WORK[I] = ZERO;
    }
  } else {
    dgemv('C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1);
  }

  dgemv('C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1);

  dgemv('N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1);
  dgemv('N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2);

  SCL.value = REALZERO;
  SSQ.value = REALZERO;
  dlassq(M1, X1, INCX1, SCL, SSQ);
  dlassq(M2, X2, INCX2, SCL, SSQ);
  NORM_NEW = SCL.value * sqrt(SSQ.value);

  // If projection is sufficiently large in norm, then stop.
  // If projection is zero, then stop.
  // Otherwise, project again.

  if (NORM_NEW >= ALPHA * NORM) {
    return;
  }

  if (NORM_NEW <= N * EPS * NORM) {
    for (IX = 1;
        INCX1 < 0 ? IX >= 1 + (M1 - 1) * INCX1 : IX <= 1 + (M1 - 1) * INCX1;
        IX += INCX1) {
      X1[IX] = ZERO;
    }
    for (IX = 1;
        INCX2 < 0 ? IX >= 1 + (M2 - 1) * INCX2 : IX <= 1 + (M2 - 1) * INCX2;
        IX += INCX2) {
      X2[IX] = ZERO;
    }
    return;
  }

  NORM = NORM_NEW;

  for (I = 1; I <= N; I++) {
    WORK[I] = ZERO;
  }

  if (M1 == 0) {
    for (I = 1; I <= N; I++) {
      WORK[I] = ZERO;
    }
  } else {
    dgemv('C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1);
  }

  dgemv('C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1);

  dgemv('N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1);
  dgemv('N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2);

  SCL.value = REALZERO;
  SSQ.value = REALZERO;
  dlassq(M1, X1, INCX1, SCL, SSQ);
  dlassq(M2, X2, INCX2, SCL, SSQ);
  NORM_NEW = SCL.value * sqrt(SSQ.value);

  // If second projection is sufficiently large in norm, then do
  // nothing more. Alternatively, if it shrunk significantly, then
  // truncate it to zero.

  if (NORM_NEW < ALPHA * NORM) {
    for (IX = 1;
        INCX1 < 0 ? IX >= 1 + (M1 - 1) * INCX1 : IX <= 1 + (M1 - 1) * INCX1;
        IX += INCX1) {
      X1[IX] = ZERO;
    }
    for (IX = 1;
        INCX2 < 0 ? IX >= 1 + (M2 - 1) * INCX2 : IX <= 1 + (M2 - 1) * INCX2;
        IX += INCX2) {
      X2[IX] = ZERO;
    }
  }
}
