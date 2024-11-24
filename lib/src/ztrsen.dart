// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/ztrexc.dart';
import 'package:lapack/src/ztrsyl.dart';

void ztrsen(
  final String JOB,
  final String COMPQ,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<Complex> W_,
  final Box<int> M,
  final Box<double> S,
  final Box<double> SEP,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final T = T_.having(ld: LDT);
  final Q = Q_.having(ld: LDQ);
  final W = W_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, WANTBH, WANTQ, WANTS, WANTSP;
  int K, KS, LWMIN = 0, N1, N2, NN;
  double RNORM;
  final ISAVE = Array<int>(3);
  final RWORK = Array<double>(1);
  final IERR = Box(0), KASE = Box(0);
  final SCALE = Box(0.0), EST = Box(0.0);

  // Decode and test the input parameters.

  WANTBH = lsame(JOB, 'B');
  WANTS = lsame(JOB, 'E') || WANTBH;
  WANTSP = lsame(JOB, 'V') || WANTBH;
  WANTQ = lsame(COMPQ, 'V');

  // Set M to the number of selected eigenvalues.

  M.value = 0;
  for (K = 1; K <= N; K++) {
    if (SELECT[K]) M.value++;
  }

  N1 = M.value;
  N2 = N - M.value;
  NN = N1 * N2;

  INFO.value = 0;
  LQUERY = (LWORK == -1);

  if (WANTSP) {
    LWMIN = max(1, 2 * NN);
  } else if (lsame(JOB, 'N')) {
    LWMIN = 1;
  } else if (lsame(JOB, 'E')) {
    LWMIN = max(1, NN);
  }

  if (!lsame(JOB, 'N') && !WANTS && !WANTSP) {
    INFO.value = -1;
  } else if (!lsame(COMPQ, 'N') && !WANTQ) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -8;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -14;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZTRSEN', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M.value == N || M.value == 0) {
    if (WANTS) S.value = ONE;
    if (WANTSP) SEP.value = zlange('1', N, N, T, LDT, RWORK);
  } else {
    // Collect the selected eigenvalues at the top left corner of T.

    KS = 0;
    for (K = 1; K <= N; K++) {
      if (SELECT[K]) {
        KS++;

        // Swap the K-th eigenvalue to position KS.

        if (K != KS) ztrexc(COMPQ, N, T, LDT, Q, LDQ, K, KS, IERR);
      }
    }

    if (WANTS) {
      // Solve the Sylvester equation for R:

      // T11*R - R*T22 = scale*T12

      zlacpy('F', N1, N2, T(1, N1 + 1), LDT, WORK.asMatrix(N1), N1);
      ztrsyl('N', 'N', -1, N1, N2, T, LDT, T(N1 + 1, N1 + 1), LDT,
          WORK.asMatrix(N1), N1, SCALE, IERR);

      // Estimate the reciprocal of the condition number of the cluster
      // of eigenvalues.

      RNORM = zlange('F', N1, N2, WORK.asMatrix(N1), N1, RWORK);
      if (RNORM == ZERO) {
        S.value = ONE;
      } else {
        S.value = SCALE.value /
            (sqrt(SCALE.value * SCALE.value / RNORM + RNORM) * sqrt(RNORM));
      }
    }

    if (WANTSP) {
      // Estimate sep(T11,T22).

      EST.value = ZERO;
      KASE.value = 0;
      while (true) {
        zlacn2(NN, WORK(NN + 1), WORK, EST, KASE, ISAVE);
        if (KASE.value == 0) break;
        if (KASE.value == 1) {
          // Solve T11*R - R*T22 = scale*X.
          ztrsyl('N', 'N', -1, N1, N2, T, LDT, T(N1 + 1, N1 + 1), LDT,
              WORK.asMatrix(N1), N1, SCALE, IERR);
        } else {
          // Solve T11**H*R - R*T22**H = scale*X.
          ztrsyl('C', 'C', -1, N1, N2, T, LDT, T(N1 + 1, N1 + 1), LDT,
              WORK.asMatrix(N1), N1, SCALE, IERR);
        }
      }

      SEP.value = SCALE.value / EST.value;
    }
  }

  // Copy reordered eigenvalues to W.
  for (K = 1; K <= N; K++) {
    W[K] = T[K][K];
  }

  WORK[1] = LWMIN.toComplex();
}
