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
import 'package:dart_lapack/src/zgemlqt.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/ztplqt.dart';
import 'package:dart_lapack/src/ztpmlqt.dart';

void zlqt05(
  final int M,
  final int N,
  final int L,
  final int NB,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 0.0;
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');
  final K = M;
  final N2 = M + N;
  final NP1 = N > 0 ? M + 1 : 1;
  final LWORK = N2 * N2 * NB;

  // Dynamically allocate all arrays

  final A = Matrix<Complex>(M, N2),
      AF = Matrix<Complex>(M, N2),
      Q = Matrix<Complex>(N2, N2),
      R = Matrix<Complex>(N2, N2),
      RWORK = Array<double>(N2),
      WORK = Array<Complex>(LWORK),
      T = Matrix<Complex>(NB, M),
      C = Matrix<Complex>(N2, M),
      CF = Matrix<Complex>(N2, M),
      D = Matrix<Complex>(M, N2),
      DF = Matrix<Complex>(M, N2);

  // Put random stuff into A

  final LDT = NB;
  zlaset('Full', M, N2, Complex.zero, Complex.zero, A, M);
  zlaset('Full', NB, M, Complex.zero, Complex.zero, T, NB);
  for (var J = 1; J <= M; J++) {
    zlarnv(2, ISEED, M - J + 1, A(J, J).asArray());
  }
  if (N > 0) {
    for (var J = 1; J <= N - L; J++) {
      zlarnv(2, ISEED, M, A(1, min(N + M, M + 1) + J - 1).asArray());
    }
  }
  if (L > 0) {
    for (var J = 1; J <= L; J++) {
      zlarnv(2, ISEED, M - J + 1,
          A(J, min(N + M, N + M - L + 1) + J - 1).asArray());
    }
  }

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N2, A, M, AF, M);

  // Factor the matrix A in the array AF.

  ztplqt(M, N, L, NB, AF, M, AF(1, NP1), M, T, LDT, WORK, INFO);

  // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

  zlaset('Full', N2, N2, Complex.zero, Complex.one, Q, N2);
  zgemlqt('L', 'N', N2, N2, K, NB, AF, M, T, LDT, Q, N2, WORK, INFO);

  // Copy L

  zlaset('Full', N2, N2, Complex.zero, Complex.zero, R, N2);
  zlacpy('Lower', M, N2, AF, M, R, N2);

  // Compute |L - A*Q*C| / |A| and store in RESULT(1)

  zgemm('N', 'C', M, N2, N2, -Complex.one, A, M, Q, N2, Complex.one, R, N2);
  final ANORM = zlange('1', M, N2, A, M, RWORK);
  var RESID = zlange('1', M, N2, R, N2, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * ANORM * max(1, N2));
  } else {
    RESULT[1] = ZERO;
  }

  // Compute |I - Q*Q'| and store in RESULT(2)

  zlaset('Full', N2, N2, Complex.zero, Complex.one, R, N2);
  zherk('U', 'N', N2, N2, -ONE, Q, N2, ONE, R, N2);
  RESID = zlansy('1', 'Upper', N2, R, N2, RWORK);
  RESULT[2] = RESID / (EPS * max(1, N2));

  // Generate random m-by-n matrix C and a copy CF

  zlaset('Full', N2, M, Complex.zero, Complex.one, C, N2);
  for (var J = 1; J <= M; J++) {
    zlarnv(2, ISEED, N2, C(1, J).asArray());
  }
  final CNORM = zlange('1', N2, M, C, N2, RWORK);
  zlacpy('Full', N2, M, C, N2, CF, N2);

  // Apply Q to C as Q*C

  ztpmlqt('L', 'N', N, M, K, L, NB, AF(1, NP1), M, T, LDT, CF, N2, CF(NP1, 1),
      N2, WORK, INFO);

  // Compute |Q*C - Q*C| / |C|

  zgemm('N', 'N', N2, M, N2, -Complex.one, Q, N2, C, N2, Complex.one, CF, N2);
  RESID = zlange('1', N2, M, CF, N2, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, N2) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  zlacpy('Full', N2, M, C, N2, CF, N2);

  // Apply Q to C as QT*C

  ztpmlqt('L', 'C', N, M, K, L, NB, AF(1, NP1), M, T, LDT, CF, N2, CF(NP1, 1),
      N2, WORK, INFO);

  // Compute |QT*C - QT*C| / |C|

  zgemm('C', 'N', N2, M, N2, -Complex.one, Q, N2, C, N2, Complex.one, CF, N2);
  RESID = zlange('1', N2, M, CF, N2, RWORK);

  if (CNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, N2) * CNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random m-by-n matrix D and a copy DF

  for (var J = 1; J <= N2; J++) {
    zlarnv(2, ISEED, M, D(1, J).asArray());
  }
  final DNORM = zlange('1', M, N2, D, M, RWORK);
  zlacpy('Full', M, N2, D, M, DF, M);

  // Apply Q to D as D*Q

  ztpmlqt('R', 'N', M, N, K, L, NB, AF(1, NP1), M, T, LDT, DF, M, DF(1, NP1), M,
      WORK, INFO);

  // Compute |D*Q - D*Q| / |D|

  zgemm('N', 'N', M, N2, N2, -Complex.one, D, M, Q, N2, Complex.one, DF, M);
  RESID = zlange('1', M, N2, DF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, N2) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  zlacpy('Full', M, N2, D, M, DF, M);

  // Apply Q to D as D*QT

  ztpmlqt('R', 'C', M, N, K, L, NB, AF(1, NP1), M, T, LDT, DF, M, DF(1, NP1), M,
      WORK, INFO);

  // Compute |D*QT - D*QT| / |D|

  zgemm('N', 'C', M, N2, N2, -Complex.one, D, M, Q, N2, Complex.one, DF, M);
  RESID = zlange('1', M, N2, DF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, N2) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
