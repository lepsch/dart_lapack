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
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/ztpmqrt.dart';
import 'package:dart_lapack/src/ztpqrt.dart';

void zqrt05(
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
  const ZERO = 0.0, ONE = 1.0;
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');
  final K = N;
  final M2 = M + N;
  final NP1 = M > 0 ? N + 1 : 1;
  final LWORK = M2 * M2 * NB;

  final A = Matrix<Complex>(M2, N),
      AF = Matrix<Complex>(M2, N),
      C = Matrix<Complex>(M2, N),
      CF = Matrix<Complex>(M2, N),
      D = Matrix<Complex>(N, M2),
      DF = Matrix<Complex>(N, M2),
      Q = Matrix<Complex>(M2, M2),
      R = Matrix<Complex>(M2, M2),
      T = Matrix<Complex>(NB, N),
      WORK = Array<Complex>(LWORK);
  final RWORK = Array<double>(M2);

  // Put random stuff into A

  final LDT = NB;
  zlaset('Full', M2, N, Complex.zero, Complex.zero, A, M2);
  zlaset('Full', NB, N, Complex.zero, Complex.zero, T, NB);
  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, J, A(1, J).asArray());
  }
  if (M > 0) {
    for (var J = 1; J <= N; J++) {
      zlarnv(2, ISEED, M - L, A(min(N + M, N + 1), J).asArray());
    }
  }
  if (L > 0) {
    for (var J = 1; J <= N; J++) {
      zlarnv(2, ISEED, min(J, L), A(min(N + M, N + M - L + 1), J).asArray());
    }
  }

  // Copy the matrix A to the array AF.

  zlacpy('Full', M2, N, A, M2, AF, M2);

  // Factor the matrix A in the array AF.

  ztpqrt(M, N, L, NB, AF, M2, AF(NP1, 1), M2, T, LDT, WORK, INFO);

  // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

  zlaset('Full', M2, M2, Complex.zero, Complex.one, Q, M2);
  zgemqrt('R', 'N', M2, M2, K, NB, AF, M2, T, LDT, Q, M2, WORK, INFO);

  // Copy R

  zlaset('Full', M2, N, Complex.zero, Complex.zero, R, M2);
  zlacpy('Upper', M2, N, AF, M2, R, M2);

  // Compute |R - Q'*A| / |A| and store in RESULT(1)

  zgemm('C', 'N', M2, N, M2, -Complex.one, Q, M2, A, M2, Complex.one, R, M2);
  final ANORM = zlange('1', M2, N, A, M2, RWORK);
  var RESID = zlange('1', M2, N, R, M2, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * ANORM * max(1, M2));
  } else {
    RESULT[1] = ZERO;
  }

  // Compute |I - Q'*Q| and store in RESULT(2)

  zlaset('Full', M2, M2, Complex.zero, Complex.one, R, M2);
  zherk('U', 'C', M2, M2, -ONE, Q, M2, ONE, R, M2);
  RESID = zlansy('1', 'Upper', M2, R, M2, RWORK);
  RESULT[2] = RESID / (EPS * max(1, M2));

  // Generate random m-by-n matrix C and a copy CF

  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M2, C(1, J).asArray());
  }
  final CNORM = zlange('1', M2, N, C, M2, RWORK);
  zlacpy('Full', M2, N, C, M2, CF, M2);

  // Apply Q to C as Q*C

  ztpmqrt('L', 'N', M, N, K, L, NB, AF(NP1, 1), M2, T, LDT, CF, M2, CF(NP1, 1),
      M2, WORK, INFO);

  // Compute |Q*C - Q*C| / |C|

  zgemm('N', 'N', M2, N, M2, -Complex.one, Q, M2, C, M2, Complex.one, CF, M2);
  RESID = zlange('1', M2, N, CF, M2, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M2) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  zlacpy('Full', M2, N, C, M2, CF, M2);

  // Apply Q to C as QT*C

  ztpmqrt('L', 'C', M, N, K, L, NB, AF(NP1, 1), M2, T, LDT, CF, M2, CF(NP1, 1),
      M2, WORK, INFO);

  // Compute |QT*C - QT*C| / |C|

  zgemm('C', 'N', M2, N, M2, -Complex.one, Q, M2, C, M2, Complex.one, CF, M2);
  RESID = zlange('1', M2, N, CF, M2, RWORK);
  if (CNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, M2) * CNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random n-by-m matrix D and a copy DF

  for (var J = 1; J <= M2; J++) {
    zlarnv(2, ISEED, N, D(1, J).asArray());
  }
  final DNORM = zlange('1', N, M2, D, N, RWORK);
  zlacpy('Full', N, M2, D, N, DF, N);

  // Apply Q to D as D*Q

  ztpmqrt('R', 'N', N, M, N, L, NB, AF(NP1, 1), M2, T, LDT, DF, N, DF(1, NP1),
      N, WORK, INFO);

  // Compute |D*Q - D*Q| / |D|

  zgemm('N', 'N', N, M2, M2, -Complex.one, D, N, Q, M2, Complex.one, DF, N);
  RESID = zlange('1', N, M2, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M2) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  zlacpy('Full', N, M2, D, N, DF, N);

  // Apply Q to D as D*QT

  ztpmqrt('R', 'C', N, M, N, L, NB, AF(NP1, 1), M2, T, LDT, DF, N, DF(1, NP1),
      N, WORK, INFO);

  // Compute |D*QT - D*QT| / |D|

  zgemm('N', 'C', N, M2, M2, -Complex.one, D, N, Q, M2, Complex.one, DF, N);
  RESID = zlange('1', N, M2, DF, N, RWORK);
  if (CNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M2) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
