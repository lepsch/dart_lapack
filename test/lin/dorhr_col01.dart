// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgemqrt.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dlatsqr.dart';
import 'package:dart_lapack/src/dorgtsqr.dart';
import 'package:dart_lapack/src/dorhr_col.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'common.dart';

void dorhr_col01(
  final int M,
  final int N,
  final int MB1,
  final int NB1,
  final int NB2,
  final Array<double> RESULT_,
) {
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 1.0;
  bool TESTZEROS;
  final WORKQUERY = Array<double>(1);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

  TESTZEROS = false;

  final EPS = dlamch('Epsilon');
  final K = min(M, N);
  final L = max(M, max(N, 1));

  // Dynamically allocate local arrays

  final A = Matrix<double>(M, N),
      AF = Matrix<double>(M, N),
      Q = Matrix<double>(L, L),
      R = Matrix<double>(M, L),
      RWORK = Array<double>(L),
      C = Matrix<double>(M, N),
      CF = Matrix<double>(M, N),
      D = Matrix<double>(N, M),
      DF = Matrix<double>(N, M);

  // Put random numbers into A and copy to AF

  for (var J = 1; J <= N; J++) {
    dlarnv(2, ISEED, M, A(1, J).asArray());
  }
  // ignore: dead_code
  if (TESTZEROS) {
    if (M >= 4) {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, M ~/ 2, A(M ~/ 4, J).asArray());
      }
    }
  }
  dlacpy('Full', M, N, A, M, AF, M);

  // Number of row blocks in DLATSQR

  final NRB = max(1, ((M - N) / (MB1 - N)).ceil()).toInt();

  final T1 = Matrix<double>(NB1, N * NRB),
      T2 = Matrix<double>(NB2, N),
      DIAG = Array<double>(N);

  // Begin determine LWORK for the array WORK and allocate memory.

  // DLATSQR requires NB1 to be bounded by N.

  final NB1_UB = min(NB1, N);

  // DGEMQRT requires NB2 to be bounded by N.

  final NB2_UB = min(NB2, N);

  dlatsqr(M, N, MB1, NB1_UB, AF, M, T1, NB1, WORKQUERY, -1, INFO);
  var LWORK = WORKQUERY[1].toInt();
  dorgtsqr(M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1, INFO);

  LWORK = max(LWORK, WORKQUERY[1].toInt());

  // In DGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
  //            or  M*NB2_UB if SIDE = 'R'.

  LWORK = max(LWORK, max(NB2_UB * N, NB2_UB * M));

  final WORK = Array<double>(LWORK);

  // End allocate memory for WORK.

  // Begin Householder reconstruction routines

  // Factor the matrix A in the array AF.

  srnamc.SRNAMT = 'DLATSQR';
  dlatsqr(M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK, INFO);

  // Copy the factor R into the array R.

  srnamc.SRNAMT = 'DLACPY';
  dlacpy('U', N, N, AF, M, R, M);

  // Reconstruct the orthogonal matrix Q.

  srnamc.SRNAMT = 'DORGTSQR';
  dorgtsqr(M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK, INFO);

  // Perform the Householder reconstruction, the result is stored
  // the arrays AF and T2.

  srnamc.SRNAMT = 'DORHR_COL';
  dorhr_col(M, N, NB2, AF, M, T2, NB2, DIAG, INFO);

  // Compute the factor R_hr corresponding to the Householder
  // reconstructed Q_hr and place it in the upper triangle of AF to
  // match the Q storage format in DGEQRT. R_hr = R_tsqr * S,
  // this means changing the sign of I-th row of the matrix R_tsqr
  // according to sign of of I-th diagonal element DIAG(I) of the
  // matrix S.

  srnamc.SRNAMT = 'DLACPY';
  dlacpy('U', N, N, R, M, AF, M);

  for (var I = 1; I <= N; I++) {
    if (DIAG[I] == -ONE) {
      dscal(N + 1 - I, -ONE, AF(I, I).asArray(), M);
    }
  }

  // End Householder reconstruction routines.

  // Generate the m-by-m matrix Q

  dlaset('Full', M, M, ZERO, ONE, Q, M);

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO);

  // Copy R

  dlaset('Full', M, N, ZERO, ZERO, R, M);

  dlacpy('Upper', M, N, AF, M, R, M);

  // TEST 1
  // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)

  dgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M);

  final ANORM = dlange('1', M, N, A, M, RWORK);
  var RESID = dlange('1', M, N, R, M, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * max(1, M) * ANORM);
  } else {
    RESULT[1] = ZERO;
  }

  // TEST 2
  // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)

  dlaset('Full', M, M, ZERO, ONE, R, M);
  dsyrk('U', 'T', M, M, -ONE, Q, M, ONE, R, M);
  RESID = dlansy('1', 'Upper', M, R, M, RWORK);
  RESULT[2] = RESID / (EPS * max(1, M));

  // Generate random m-by-n matrix C

  for (var J = 1; J <= N; J++) {
    dlarnv(2, ISEED, M, C(1, J).asArray());
  }
  final CNORM = dlange('1', M, N, C, M, RWORK);
  dlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as Q*C = CF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO);

  // TEST 3
  // Compute |CF - Q*C| / ( eps *  m * |C| )

  dgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M);
  RESID = dlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  dlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as (Q**T)*C = CF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('L', 'T', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO);

  // TEST 4
  // Compute |CF - (Q**T)*C| / ( eps * m * |C|)

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

  // Apply Q to D as D*Q = DF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO);

  // TEST 5
  // Compute |DF - D*Q| / ( eps * m * |D| )

  dgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N);
  RESID = dlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  dlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*QT = DF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('R', 'T', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO);

  // TEST 6
  // Compute |DF - D*(Q**T)| / ( eps * m * |D| )

  dgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N);
  RESID = dlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
