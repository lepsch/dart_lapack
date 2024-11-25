// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgemqrt.dart';
import 'package:dart_lapack/src/dgeqrt.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dqrt04(
  final int M,
  final int N,
  final int NB,
  final Array<double> RESULT_,
) {
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 1.0;
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');
  final K = min(M, N);
  final L = max(M, N);
  final LWORK = max(2, L) * max(2, L).toInt() * NB;

  // Dynamically allocate local arrays

  final A = Matrix<double>(M, N),
      AF = Matrix<double>(M, N),
      Q = Matrix<double>(M, M),
      R = Matrix<double>(M, L),
      RWORK = Array<double>(L),
      WORK = Array<double>(LWORK),
      T = Matrix<double>(NB, N),
      C = Matrix<double>(M, N),
      CF = Matrix<double>(M, N),
      D = Matrix<double>(N, M),
      DF = Matrix<double>(N, M);

  // Put random numbers into A and copy to AF

  final LDT = NB;
  for (var J = 1; J <= N; J++) {
    dlarnv(2, ISEED, M, A(1, J).asArray());
  }
  dlacpy('Full', M, N, A, M, AF, M);

  // Factor the matrix A in the array AF.

  dgeqrt(M, N, NB, AF, M, T, LDT, WORK, INFO);

  // Generate the m-by-m matrix Q

  dlaset('Full', M, M, ZERO, ONE, Q, M);
  dgemqrt('R', 'N', M, M, K, NB, AF, M, T, LDT, Q, M, WORK, INFO);

  // Copy R

  dlaset('Full', M, N, ZERO, ZERO, R, M);
  dlacpy('Upper', M, N, AF, M, R, M);

  // Compute |R - Q'*A| / |A| and store in RESULT(1)

  dgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M);
  final ANORM = dlange('1', M, N, A, M, RWORK);
  var RESID = dlange('1', M, N, R, M, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * max(1, M) * ANORM);
  } else {
    RESULT[1] = ZERO;
  }

  // Compute |I - Q'*Q| and store in RESULT(2)

  dlaset('Full', M, M, ZERO, ONE, R, M);
  dsyrk('U', 'C', M, M, -ONE, Q, M, ONE, R, M);
  RESID = dlansy('1', 'Upper', M, R, M, RWORK);
  RESULT[2] = RESID / (EPS * max(1, M));

  // Generate random m-by-n matrix C and a copy CF

  for (var J = 1; J <= N; J++) {
    dlarnv(2, ISEED, M, C(1, J).asArray());
  }
  final CNORM = dlange('1', M, N, C, M, RWORK);
  dlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as Q*C

  dgemqrt('L', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

  // Compute |Q*C - Q*C| / |C|

  dgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M);
  RESID = dlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  dlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as QT*C

  dgemqrt('L', 'T', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

  // Compute |QT*C - QT*C| / |C|

  dgemm('T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M);
  RESID = dlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random n-by-m matrix D and a copy DF

  for (var J = 1; J <= M; J++) {
    dlarnv(2, ISEED, N, D(1, J).asArray());
  }
  final DNORM = dlange('1', N, M, D, N, RWORK);
  dlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*Q

  dgemqrt('R', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

  // Compute |D*Q - D*Q| / |D|

  dgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N);
  RESID = dlange('1', N, M, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  dlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*QT

  dgemqrt('R', 'T', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

  // Compute |D*QT - D*QT| / |D|

  dgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N);
  RESID = dlange('1', N, M, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
