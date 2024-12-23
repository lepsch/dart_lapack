// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zunt01(
  final String ROWCOL,
  final int M,
  final int N,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  String TRANSU;
  int I, J, K, LDWORK, MNMIN;
  double EPS;
  Complex TMP;

  RESID.value = ZERO;

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  EPS = dlamch('Precision');
  if (M < N || (M == N && lsame(ROWCOL, 'R'))) {
    TRANSU = 'N';
    K = N;
  } else {
    TRANSU = 'C';
    K = M;
  }
  MNMIN = min(M, N);

  if ((MNMIN + 1) * MNMIN <= LWORK) {
    LDWORK = MNMIN;
  } else {
    LDWORK = 0;
  }
  if (LDWORK > 0) {
    // Compute I - U*U' or I - U'*U.

    zlaset('Upper', MNMIN, MNMIN, Complex.zero, Complex.one,
        WORK.asMatrix(LDWORK), LDWORK);
    zherk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK.asMatrix(LDWORK),
        LDWORK);

    // Compute norm( I - U*U' ) / ( K * EPS ) .

    RESID.value =
        zlansy('1', 'Upper', MNMIN, WORK.asMatrix(LDWORK), LDWORK, RWORK);
    RESID.value = (RESID.value / K) / EPS;
  } else if (TRANSU == 'C') {
    // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        if (I != J) {
          TMP = Complex.zero;
        } else {
          TMP = Complex.one;
        }
        TMP -= zdotc(M, U(1, I).asArray(), 1, U(1, J).asArray(), 1);
        RESID.value = max(RESID.value, TMP.cabs1());
      }
    }
    RESID.value = (RESID.value / M) / EPS;
  } else {
    // Find the maximum element in abs( I - U*U' ) / ( n * EPS )

    for (J = 1; J <= M; J++) {
      for (I = 1; I <= J; I++) {
        if (I != J) {
          TMP = Complex.zero;
        } else {
          TMP = Complex.one;
        }
        TMP -= zdotc(N, U(J, 1).asArray(), LDU, U(I, 1).asArray(), LDU);
        RESID.value = max(RESID.value, TMP.cabs1());
      }
    }
    RESID.value = (RESID.value / N) / EPS;
  }
}
