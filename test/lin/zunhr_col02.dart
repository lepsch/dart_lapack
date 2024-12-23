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
import 'package:dart_lapack/src/zgetsqrhrt.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';

import 'common.dart';

void zunhr_col02(
  final int M,
  final int N,
  final int MB1,
  final int NB1,
  final int NB2,
  final Array<double> RESULT_,
) {
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 1.0;
  final WORKQUERY = Array<Complex>(1);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

  final TESTZEROS = false;

  final EPS = dlamch('Epsilon');
  final K = min(M, N);
  final L = max(M, max(N, 1));

  // Dynamically allocate local arrays

  final A = Matrix<Complex>(M, N),
      AF = Matrix<Complex>(M, N),
      Q = Matrix<Complex>(L, L),
      R = Matrix<Complex>(M, L),
      RWORK = Array<double>(L),
      C = Matrix<Complex>(M, N),
      CF = Matrix<Complex>(M, N),
      D = Matrix<Complex>(N, M),
      DF = Matrix<Complex>(N, M);

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

  // Number of row blocks in ZLATSQR

  // final NRB = max(1, ((M - N) / (MB1 - N)).ceil());

  // final T1 = Matrix<Complex>(NB1, N * NRB);
  final T2 = Matrix<Complex>(NB2, N);
  // final DIAG = Array<Complex>(N);

  // Begin determine LWORK for the array WORK and allocate memory.

  // ZGEMQRT requires NB2 to be bounded by N.

  final NB2_UB = min(NB2, N);

  zgetsqrhrt(M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORKQUERY, -1, INFO);

  var LWORK = WORKQUERY[1].toInt();

  // In ZGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
  //            or  M*NB2_UB if SIDE = 'R'.

  LWORK = max(LWORK, max(NB2_UB * N, NB2_UB * M));

  final WORK = Array<Complex>(LWORK);

  // End allocate memory for WORK.

  // Begin Householder reconstruction routines

  // Factor the matrix A in the array AF.

  srnamc.SRNAMT = 'ZGETSQRHRT';
  zgetsqrhrt(M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORK, LWORK, INFO);

  // End Householder reconstruction routines.

  // Generate the m-by-m matrix Q

  zlaset('Full', M, M, Complex.zero, Complex.one, Q, M);

  srnamc.SRNAMT = 'ZGEMQRT';
  zgemqrt('L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO);

  // Copy R

  zlaset('Full', M, N, Complex.zero, Complex.zero, R, M);

  zlacpy('Upper', M, N, AF, M, R, M);

  // TEST 1
  // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)

  zgemm('C', 'N', M, N, M, -Complex.one, Q, M, A, M, Complex.one, R, M);

  final ANORM = zlange('1', M, N, A, M, RWORK);
  var RESID = zlange('1', M, N, R, M, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * max(1, M) * ANORM);
  } else {
    RESULT[1] = ZERO;
  }

  // TEST 2
  // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)

  zlaset('Full', M, M, Complex.zero, Complex.one, R, M);
  zherk('U', 'C', M, M, -ONE, Q, M, ONE, R, M);
  RESID = zlansy('1', 'Upper', M, R, M, RWORK);
  RESULT[2] = RESID / (EPS * max(1, M));

  // Generate random m-by-n matrix C

  for (var J = 1; J <= N; J++) {
    zlarnv(2, ISEED, M, C(1, J).asArray());
  }
  final CNORM = zlange('1', M, N, C, M, RWORK);
  zlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as Q*C = CF

  srnamc.SRNAMT = 'ZGEMQRT';
  zgemqrt('L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO);

  // TEST 3
  // Compute |CF - Q*C| / ( eps *  m * |C| )

  zgemm('N', 'N', M, N, M, -Complex.one, Q, M, C, M, Complex.one, CF, M);
  RESID = zlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  zlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as (Q**T)*C = CF

  srnamc.SRNAMT = 'ZGEMQRT';
  zgemqrt('L', 'C', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO);

  // TEST 4
  // Compute |CF - (Q**T)*C| / ( eps * m * |C|)

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

  // Apply Q to D as D*Q = DF

  srnamc.SRNAMT = 'ZGEMQRT';
  zgemqrt('R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO);

  // TEST 5
  // Compute |DF - D*Q| / ( eps * m * |D| )

  zgemm('N', 'N', N, M, M, -Complex.one, D, N, Q, M, Complex.one, DF, N);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  zlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*QT = DF

  srnamc.SRNAMT = 'ZGEMQRT';
  zgemqrt('R', 'C', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO);

  // TEST 6
  // Compute |DF - D*(Q**T)| / ( eps * m * |D| )

  zgemm('N', 'C', N, M, M, -Complex.one, D, N, Q, M, Complex.one, DF, N);
  RESID = zlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
