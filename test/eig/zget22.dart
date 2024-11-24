// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zget22(
  final String TRANSA,
  final String TRANSE,
  final String TRANSW,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> E_,
  final int LDE,
  final Array<Complex> W_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having(ld: LDE);
  final W = W_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having(length: 2);

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  String NORMA, NORME;
  int ITRNSE, ITRNSW, J, JCOL, JOFF, JROW, JVEC;
  double ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL;
  Complex WTEMP;

  // Initialize RESULT (in case N=0)

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');

  ITRNSE = 0;
  ITRNSW = 0;
  NORMA = 'O';
  NORME = 'O';

  if (lsame(TRANSA, 'T') || lsame(TRANSA, 'C')) {
    NORMA = 'I';
  }

  if (lsame(TRANSE, 'T')) {
    ITRNSE = 1;
    NORME = 'I';
  } else if (lsame(TRANSE, 'C')) {
    ITRNSE = 2;
    NORME = 'I';
  }

  if (lsame(TRANSW, 'C')) {
    ITRNSW = 1;
  }

  // Normalization of E:

  ENRMIN = ONE / ULP;
  ENRMAX = ZERO;
  if (ITRNSE == 0) {
    for (JVEC = 1; JVEC <= N; JVEC++) {
      TEMP1 = ZERO;
      for (J = 1; J <= N; J++) {
        TEMP1 = max(TEMP1, E[J][JVEC].real.abs() + E[J][JVEC].imaginary.abs());
      }
      ENRMIN = min(ENRMIN, TEMP1);
      ENRMAX = max(ENRMAX, TEMP1);
    }
  } else {
    for (JVEC = 1; JVEC <= N; JVEC++) {
      RWORK[JVEC] = ZERO;
    }

    for (J = 1; J <= N; J++) {
      for (JVEC = 1; JVEC <= N; JVEC++) {
        RWORK[JVEC] = max(
            RWORK[JVEC], E[JVEC][J].real.abs() + E[JVEC][J].imaginary.abs());
      }
    }

    for (JVEC = 1; JVEC <= N; JVEC++) {
      ENRMIN = min(ENRMIN, RWORK[JVEC]);
      ENRMAX = max(ENRMAX, RWORK[JVEC]);
    }
  }

  // Norm of A:

  ANORM = max(zlange(NORMA, N, N, A, LDA, RWORK), UNFL);

  // Norm of E:

  ENORM = max(zlange(NORME, N, N, E, LDE, RWORK), ULP);

  // Norm of error:

  // Error =  AE - EW

  zlaset('Full', N, N, Complex.zero, Complex.zero, WORK.asMatrix(), N);

  JOFF = 0;
  for (JCOL = 1; JCOL <= N; JCOL++) {
    if (ITRNSW == 0) {
      WTEMP = W[JCOL];
    } else {
      WTEMP = W[JCOL].conjugate();
    }

    if (ITRNSE == 0) {
      for (JROW = 1; JROW <= N; JROW++) {
        WORK[JOFF + JROW] = E[JROW][JCOL] * WTEMP;
      }
    } else if (ITRNSE == 1) {
      for (JROW = 1; JROW <= N; JROW++) {
        WORK[JOFF + JROW] = E[JCOL][JROW] * WTEMP;
      }
    } else {
      for (JROW = 1; JROW <= N; JROW++) {
        WORK[JOFF + JROW] = E[JCOL][JROW].conjugate() * WTEMP;
      }
    }
    JOFF += N;
  }

  zgemm(TRANSA, TRANSE, N, N, N, Complex.one, A, LDA, E, LDE, -Complex.one,
      WORK.asMatrix(), N);

  ERRNRM = zlange('One', N, N, WORK.asMatrix(), N, RWORK) / ENORM;

  // Compute RESULT(1) (avoiding under/overflow)

  if (ANORM > ERRNRM) {
    RESULT[1] = (ERRNRM / ANORM) / ULP;
  } else {
    if (ANORM < ONE) {
      RESULT[1] = ONE / ULP;
    } else {
      RESULT[1] = min(ERRNRM / ANORM, ONE) / ULP;
    }
  }

  // Compute RESULT(2) : the normalization error in E.

  RESULT[2] = max((ENRMAX - ONE).abs(), (ENRMIN - ONE).abs()) / (N * ULP);
}
