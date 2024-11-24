// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgemqrt.dart';
import 'package:dart_lapack/src/zgeqrt.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zqrt04(
  final int M,
  final int N,
  final int NB,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 1.0;
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final EPS = dlamch('Epsilon');
  final K = min(M, N);
  final L = max(M, N);
  final LWORK = max(2, L) * max(2, L).toInt() * NB;
  final A = Matrix<Complex>(M, N),
      AF = Matrix<Complex>(M, N),
      Q = Matrix<Complex>(M, M),
      R = Matrix<Complex>(M, L),
      RWORK = Array<double>(L),
      WORK = Array<Complex>(LWORK),
      T = Matrix<Complex>(NB, N),
      C = Matrix<Complex>(M, N),
      CF = Matrix<Complex>(M, N),
      D = Matrix<Complex>(N, M),
      DF = Matrix<Complex>(N, M);
  final INFO = Box(0);

  // Put random numbers into A and copy to AF

  final LDT = NB;
  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M, A(1, J).asArray());
  }
  zlacpy('Full', M, N, A, M, AF, M);

  // Factor the matrix A in the array AF.

  zgeqrt(M, N, NB, AF, M, T, LDT, WORK, INFO);

  // Generate the m-by-m matrix Q

  zlaset('Full', M, M, Complex.zero, Complex.one, Q, M);
  zgemqrt('R', 'N', M, M, K, NB, AF, M, T, LDT, Q, M, WORK, INFO);

  // Copy R

  zlaset('Full', M, N, Complex.zero, Complex.zero, R, M);
  zlacpy('Upper', M, N, AF, M, R, M);

  // Compute |R - Q'*A| / |A| and store in RESULT(1)

  zgemm('C', 'N', M, N, M, -Complex.one, Q, M, A, M, Complex.one, R, M);
  final ANORM = zlange('1', M, N, A, M, RWORK);
  var RESID = zlange('1', M, N, R, M, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * max(1, M) * ANORM);
  } else {
    RESULT[1] = ZERO;
  }

  // Compute |I - Q'*Q| and store in RESULT(2)

  zlaset('Full', M, M, Complex.zero, Complex.one, R, M);
  zherk('U', 'C', M, M, -ONE, Q, M, ONE, R, M);
  RESID = zlansy('1', 'Upper', M, R, M, RWORK);
  RESULT[2] = RESID / (EPS * max(1, M));

  // Generate random m-by-n matrix C and a copy CF

  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M, C(1, J).asArray());
  }
  final CNORM = zlange('1', M, N, C, M, RWORK);
  zlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as Q*C

  zgemqrt('L', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

  // Compute |Q*C - Q*C| / |C|

  zgemm('N', 'N', M, N, M, -Complex.one, Q, M, C, M, Complex.one, CF, M);
  RESID = zlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  zlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as QT*C

  zgemqrt('L', 'C', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

  // Compute |QT*C - QT*C| / |C|

  zgemm('C', 'N', M, N, M, -Complex.one, Q, M, C, M, Complex.one, CF, M);
  RESID = zlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random n-by-m matrix D and a copy DF

  for (var J = 1; J <= M; J++) {
    zlarnv(2, ISEED, N, D(1, J).asArray());
  }
  final DNORM = zlange('1', N, M, D, N, RWORK);
  zlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*Q

  zgemqrt('R', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

  // Compute |D*Q - D*Q| / |D|

  zgemm('N', 'N', N, M, M, -Complex.one, D, N, Q, M, Complex.one, DF, N);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  zlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*QT

  zgemqrt('R', 'C', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

  // Compute |D*QT - D*QT| / |D|

  zgemm('N', 'C', N, M, M, -Complex.one, D, N, Q, M, Complex.one, DF, N);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
