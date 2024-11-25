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
import 'package:dart_lapack/src/zgelqt.dart';
import 'package:dart_lapack/src/zgemlqt.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zlqt04(
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
  final LL = max(M, N);
  final LWORK = max(2, LL) * max(2, LL).toInt() * NB;

  final A = Matrix<Complex>(M, N),
      AF = Matrix<Complex>(M, N),
      Q = Matrix<Complex>(N, N),
      L = Matrix<Complex>(LL, N),
      RWORK = Array<double>(LL),
      WORK = Array<Complex>(LWORK),
      T = Matrix<Complex>(NB, N),
      C = Matrix<Complex>(M, N),
      CF = Matrix<Complex>(M, N),
      D = Matrix<Complex>(N, M),
      DF = Matrix<Complex>(N, M);

  // Put random numbers into A and copy to AF

  final LDT = NB;
  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M, A(1, J).asArray());
  }
  zlacpy('Full', M, N, A, M, AF, M);

  // Factor the matrix A in the array AF.

  zgelqt(M, N, NB, AF, M, T, LDT, WORK, INFO);

  // Generate the n-by-n matrix Q

  zlaset('Full', N, N, Complex.zero, Complex.one, Q, N);
  zgemlqt('R', 'N', N, N, K, NB, AF, M, T, LDT, Q, N, WORK, INFO);

  // Copy L

  zlaset('Full', LL, N, Complex.zero, Complex.zero, L, LL);
  zlacpy('Lower', M, N, AF, M, L, LL);

  // Compute |L - A*Q'| / |A| and store in RESULT(1)

  zgemm('N', 'C', M, N, N, -Complex.one, A, M, Q, N, Complex.one, L, LL);
  final ANORM = zlange('1', M, N, A, M, RWORK);
  var RESID = zlange('1', M, N, L, LL, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * max(1, M) * ANORM);
  } else {
    RESULT[1] = ZERO;
  }

  // Compute |I - Q'*Q| and store in RESULT(2)

  zlaset('Full', N, N, Complex.zero, Complex.one, L, LL);
  zherk('U', 'C', N, N, -ONE, Q, N, ONE, L, LL);
  RESID = zlansy('1', 'Upper', N, L, LL, RWORK);
  RESULT[2] = RESID / (EPS * max(1, N));

  // Generate random m-by-n matrix C and a copy CF

  for (var J = 1; J <= M; J++) {
    zlarnv(2, ISEED, N, D(1, J).asArray());
  }
  final DNORM = zlange('1', N, M, D, N, RWORK);
  zlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to C as Q*C

  zgemlqt('L', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

  // Compute |Q*D - Q*D| / |D|

  zgemm('N', 'N', N, M, N, -Complex.one, Q, N, D, N, Complex.one, DF, N);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy D into DF again

  zlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as QT*D

  zgemlqt('L', 'C', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

  // Compute |QT*D - QT*D| / |D|

  zgemm('C', 'N', N, M, N, -Complex.one, Q, N, D, N, Complex.one, DF, N);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random n-by-m matrix D and a copy DF

  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M, C(1, J).asArray());
  }
  final CNORM = zlange('1', M, N, C, M, RWORK);
  zlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as C*Q

  zgemlqt('R', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

  // Compute |C*Q - C*Q| / |C|

  zgemm('N', 'N', M, N, N, -Complex.one, C, M, Q, N, Complex.one, CF, M);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy C into CF again

  zlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to D as D*QT

  zgemlqt('R', 'C', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

  // Compute |C*QT - C*QT| / |C|

  zgemm('N', 'C', M, N, N, -Complex.one, C, M, Q, N, Complex.one, CF, M);
  RESID = zlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
