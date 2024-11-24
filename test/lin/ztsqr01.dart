// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgelq.dart';
import 'package:dart_lapack/src/zgemlq.dart';
import 'package:dart_lapack/src/zgemqr.dart';
import 'package:dart_lapack/src/zgeqr.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';

import 'common.dart';

void ztsqr01(
  final String TSSW,
  final int M,
  final int N,
  final int MB,
  final int NB,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 1.0;
  final TQUERY = Array<Complex>(5), WORKQUERY = Array<Complex>(1);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  // TEST TALL SKINNY OR SHORT WIDE

  final TS = lsame(TSSW, 'TS');

  // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

  var TESTZEROS = false;

  final EPS = dlamch('Epsilon');
  final K = min(M, N);
  final L = max(M, max(N, 1));
  // final MNB = max(MB, NB);
  // final LWORK = max(3, L) * MNB;

  final A = Matrix<Complex>(M, N),
      AF = Matrix<Complex>(M, N),
      Q = Matrix<Complex>(L, L),
      R = Matrix<Complex>(M, L),
      RWORK = Array<double>(L),
      C = Matrix<Complex>(M, N),
      CF = Matrix<Complex>(M, N),
      D = Matrix<Complex>(N, M),
      DF = Matrix<Complex>(N, M),
      LQ = Matrix<Complex>(L, N);

  // Put random numbers into A and copy to AF

  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M, A(1, J).asArray());
  }
  // ignore: dead_code
  if (TESTZEROS) {
    if (M >= 4) {
      for (var J = 1; J <= N; J++) {
        zlarnv(2, ISEED, M ~/ 2, A(M ~/ 4, J).asArray());
      }
    }
  }
  zlacpy('Full', M, N, A, M, AF, M);

  if (TS) {
    // Factor the matrix A in the array AF.

    zgeqr(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO);
    final TSIZE = TQUERY[1].toInt();
    var LWORK = WORKQUERY[1].toInt();
    zgemqr('L', 'N', M, M, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemqr('L', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemqr('L', 'C', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemqr('R', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemqr('R', 'C', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    final T = Array<Complex>(TSIZE);
    final WORK = Array<Complex>(LWORK);
    srnamc.SRNAMT = 'ZGEQR';
    zgeqr(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO);

    // Generate the m-by-m matrix Q

    zlaset('Full', M, M, Complex.zero, Complex.one, Q, M);
    srnamc.SRNAMT = 'ZGEMQR';
    zgemqr('L', 'N', M, M, K, AF, M, T, TSIZE, Q, M, WORK, LWORK, INFO);

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

    srnamc.SRNAMT = 'ZGEMQR';
    zgemqr('L', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

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

    srnamc.SRNAMT = 'ZGEMQR';
    zgemqr('L', 'C', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

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

    srnamc.SRNAMT = 'ZGEMQR';
    zgemqr('R', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

    // Compute |D*Q - D*Q| / |D|

    zgemm('N', 'N', N, M, M, -Complex.one, D, N, Q, M, Complex.one, DF, N);
    RESID = zlange('1', N, M, DF, N, RWORK);
    if (DNORM > ZERO) {
      RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
    } else {
      RESULT[5] = ZERO;
    }

    // Copy D into DF again

    zlacpy('Full', N, M, D, N, DF, N);

    // Apply Q to D as D*QT

    zgemqr('R', 'C', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

    // Compute |D*QT - D*QT| / |D|

    zgemm('N', 'C', N, M, M, -Complex.one, D, N, Q, M, Complex.one, DF, N);
    RESID = zlange('1', N, M, DF, N, RWORK);
    if (CNORM > ZERO) {
      RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
    } else {
      RESULT[6] = ZERO;
    }

    // Short and wide
  } else {
    zgelq(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO);
    final TSIZE = TQUERY[1].toInt();
    var LWORK = WORKQUERY[1].toInt();
    zgemlq('R', 'N', N, N, K, AF, M, TQUERY, TSIZE, Q, N, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemlq('L', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemlq('L', 'C', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemlq('R', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    zgemlq('R', 'C', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
    LWORK = max(LWORK, WORKQUERY[1].toInt());
    final T = Array<Complex>(TSIZE);
    final WORK = Array<Complex>(LWORK);
    srnamc.SRNAMT = 'ZGELQ';
    zgelq(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO);

    // Generate the n-by-n matrix Q

    zlaset('Full', N, N, Complex.zero, Complex.one, Q, N);
    srnamc.SRNAMT = 'ZGEMLQ';
    zgemlq('R', 'N', N, N, K, AF, M, T, TSIZE, Q, N, WORK, LWORK, INFO);

    // Copy R

    zlaset('Full', M, N, Complex.zero, Complex.zero, LQ, L);
    zlacpy('Lower', M, N, AF, M, LQ, L);

    // Compute |L - A*Q'| / |A| and store in RESULT(1)

    zgemm('N', 'C', M, N, N, -Complex.one, A, M, Q, N, Complex.one, LQ, L);
    final ANORM = zlange('1', M, N, A, M, RWORK);
    var RESID = zlange('1', M, N, LQ, L, RWORK);
    if (ANORM > ZERO) {
      RESULT[1] = RESID / (EPS * max(1, N) * ANORM);
    } else {
      RESULT[1] = ZERO;
    }

    // Compute |I - Q'*Q| and store in RESULT(2)

    zlaset('Full', N, N, Complex.zero, Complex.one, LQ, L);
    zherk('U', 'C', N, N, -ONE, Q, N, ONE, LQ, L);
    RESID = zlansy('1', 'Upper', N, LQ, L, RWORK);
    RESULT[2] = RESID / (EPS * max(1, N));

    // Generate random m-by-n matrix C and a copy CF

    for (var J = 1; J <= M; J++) {
      zlarnv(2, ISEED, N, D(1, J).asArray());
    }
    final DNORM = zlange('1', N, M, D, N, RWORK);
    zlacpy('Full', N, M, D, N, DF, N);

    // Apply Q to C as Q*C

    zgemlq('L', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

    // Compute |Q*D - Q*D| / |D|

    zgemm('N', 'N', N, M, N, -Complex.one, Q, N, D, N, Complex.one, DF, N);
    RESID = zlange('1', N, M, DF, N, RWORK);
    if (DNORM > ZERO) {
      RESULT[3] = RESID / (EPS * max(1, N) * DNORM);
    } else {
      RESULT[3] = ZERO;
    }

    // Copy D into DF again

    zlacpy('Full', N, M, D, N, DF, N);

    // Apply Q to D as QT*D

    zgemlq('L', 'C', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

    // Compute |QT*D - QT*D| / |D|

    zgemm('C', 'N', N, M, N, -Complex.one, Q, N, D, N, Complex.one, DF, N);
    RESID = zlange('1', N, M, DF, N, RWORK);
    if (DNORM > ZERO) {
      RESULT[4] = RESID / (EPS * max(1, N) * DNORM);
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

    zgemlq('R', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

    // Compute |C*Q - C*Q| / |C|

    zgemm('N', 'N', M, N, N, -Complex.one, C, M, Q, N, Complex.one, CF, M);
    RESID = zlange('1', N, M, DF, N, RWORK);
    if (CNORM > ZERO) {
      RESULT[5] = RESID / (EPS * max(1, N) * CNORM);
    } else {
      RESULT[5] = ZERO;
    }

    // Copy C into CF again

    zlacpy('Full', M, N, C, M, CF, M);

    // Apply Q to D as D*QT

    zgemlq('R', 'C', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

    // Compute |C*QT - C*QT| / |C|

    zgemm('N', 'C', M, N, N, -Complex.one, C, M, Q, N, Complex.one, CF, M);
    RESID = zlange('1', M, N, CF, M, RWORK);
    if (CNORM > ZERO) {
      RESULT[6] = RESID / (EPS * max(1, N) * CNORM);
    } else {
      RESULT[6] = ZERO;
    }
  }
}
