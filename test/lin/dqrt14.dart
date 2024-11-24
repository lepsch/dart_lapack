// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelq2.dart';
import 'package:lapack/src/dgeqr2.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double dqrt14(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Array<double> WORK_,
  final int LWORK,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having(length: LWORK);

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  final int LDWORK;
  final bool TPSD;
  if (lsame(TRANS, 'N')) {
    LDWORK = M + NRHS;
    TPSD = false;
    if (LWORK < (M + NRHS) * (N + 2)) {
      xerbla('DQRT14', 10);
      return ZERO;
    } else if (N <= 0 || NRHS <= 0) {
      return ZERO;
    }
  } else if (lsame(TRANS, 'T')) {
    LDWORK = M;
    TPSD = true;
    if (LWORK < (N + NRHS) * (M + 2)) {
      xerbla('DQRT14', 10);
      return ZERO;
    } else if (M <= 0 || NRHS <= 0) {
      return ZERO;
    }
  } else {
    xerbla('DQRT14', 1);
    return ZERO;
  }

  // Copy and scale A

  dlacpy('All', M, N, A, LDA, WORK.asMatrix(), LDWORK);
  final ANRM = dlange('M', M, N, WORK.asMatrix(), LDWORK, RWORK);
  if (ANRM != ZERO) {
    dlascl('G', 0, 0, ANRM, ONE, M, N, WORK.asMatrix(), LDWORK, INFO);
  }

  // Copy X or X' into the right place and scale it

  double ERR;
  if (TPSD) {
    // Copy X into columns n+1:n+nrhs of work

    dlacpy('All', M, NRHS, X, LDX, WORK(N * LDWORK + 1).asMatrix(), LDWORK);
    final XNRM =
        dlange('M', M, NRHS, WORK(N * LDWORK + 1).asMatrix(), LDWORK, RWORK);
    if (XNRM != ZERO) {
      dlascl('G', 0, 0, XNRM, ONE, M, NRHS, WORK(N * LDWORK + 1).asMatrix(),
          LDWORK, INFO);
    }

    // Compute QR factorization of X

    dgeqr2(M, N + NRHS, WORK.asMatrix(), LDWORK, WORK(LDWORK * (N + NRHS) + 1),
        WORK(LDWORK * (N + NRHS) + min(M, N + NRHS).toInt() + 1), INFO);

    // Compute largest entry in upper triangle of
    // work(n+1:m,n+1:n+nrhs)

    ERR = ZERO;
    for (var J = N + 1; J <= N + NRHS; J++) {
      for (var I = N + 1; I <= min(M, J); I++) {
        ERR = max(ERR, WORK[I + (J - 1) * M].abs());
      }
    }
  } else {
    // Copy X' into rows m+1:m+nrhs of work

    for (var I = 1; I <= N; I++) {
      for (var J = 1; J <= NRHS; J++) {
        WORK[M + J + (I - 1) * LDWORK] = X[I][J];
      }
    }

    final XNRM = dlange('M', NRHS, N, WORK(M + 1).asMatrix(), LDWORK, RWORK);
    if (XNRM != ZERO) {
      dlascl(
          'G', 0, 0, XNRM, ONE, NRHS, N, WORK(M + 1).asMatrix(), LDWORK, INFO);
    }

    // Compute LQ factorization of work

    dgelq2(LDWORK, N, WORK.asMatrix(), LDWORK, WORK(LDWORK * N + 1),
        WORK(LDWORK * (N + 1) + 1), INFO);

    // Compute largest entry in lower triangle in
    // work(m+1:m+nrhs,m+1:n)

    ERR = ZERO;
    for (var J = M + 1; J <= N; J++) {
      for (var I = J; I <= LDWORK; I++) {
        ERR = max(ERR, WORK[I + (J - 1) * LDWORK].abs());
      }
    }
  }

  return ERR / (max(M, max(N, NRHS)) * dlamch('Epsilon'));
}
