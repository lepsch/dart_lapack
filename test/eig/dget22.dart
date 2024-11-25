// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dget22(
  final String TRANSA,
  final String TRANSE,
  final String TRANSW,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> E_,
  final int LDE,
  final Array<double> WR_,
  final Array<double> WI_,
  final Array<double> WORK_,
  final Array<double> RESULT,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having(ld: LDE);
  final WR = WR_.having();
  final WI = WI_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  String NORMA, NORME;
  int IECOL, IEROW, INCE, IPAIR, ITRNSE, J, JCOL, JVEC;
  double ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL;
  final WMAT = Matrix<double>(2, 2);

  // Initialize RESULT (in case N=0)

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');

  ITRNSE = 0;
  INCE = 1;
  NORMA = 'O';
  NORME = 'O';

  if (lsame(TRANSA, 'T') || lsame(TRANSA, 'C')) {
    NORMA = 'I';
  }
  if (lsame(TRANSE, 'T') || lsame(TRANSE, 'C')) {
    NORME = 'I';
    ITRNSE = 1;
    INCE = LDE;
  }

  // Check normalization of E

  ENRMIN = ONE / ULP;
  ENRMAX = ZERO;
  if (ITRNSE == 0) {
    // Eigenvectors are column vectors.

    IPAIR = 0;
    for (JVEC = 1; JVEC <= N; JVEC++) {
      TEMP1 = ZERO;
      if (IPAIR == 0 && JVEC < N && WI[JVEC] != ZERO) IPAIR = 1;
      if (IPAIR == 1) {
        // Complex eigenvector

        for (J = 1; J <= N; J++) {
          TEMP1 = max(TEMP1, E[J][JVEC].abs() + E[J][JVEC + 1].abs());
        }
        ENRMIN = min(ENRMIN, TEMP1);
        ENRMAX = max(ENRMAX, TEMP1);
        IPAIR = 2;
      } else if (IPAIR == 2) {
        IPAIR = 0;
      } else {
        // Real eigenvector

        for (J = 1; J <= N; J++) {
          TEMP1 = max(TEMP1, E[J][JVEC].abs());
        }
        ENRMIN = min(ENRMIN, TEMP1);
        ENRMAX = max(ENRMAX, TEMP1);
        IPAIR = 0;
      }
    }
  } else {
    // Eigenvectors are row vectors.

    for (JVEC = 1; JVEC <= N; JVEC++) {
      WORK[JVEC] = ZERO;
    }

    for (J = 1; J <= N; J++) {
      IPAIR = 0;
      for (JVEC = 1; JVEC <= N; JVEC++) {
        if (IPAIR == 0 && JVEC < N && WI[JVEC] != ZERO) IPAIR = 1;
        if (IPAIR == 1) {
          WORK[JVEC] = max(WORK[JVEC], E[J][JVEC].abs() + E[J][JVEC + 1].abs());
          WORK[JVEC + 1] = WORK[JVEC];
        } else if (IPAIR == 2) {
          IPAIR = 0;
        } else {
          WORK[JVEC] = max(WORK[JVEC], E[J][JVEC].abs());
          IPAIR = 0;
        }
      }
    }

    for (JVEC = 1; JVEC <= N; JVEC++) {
      ENRMIN = min(ENRMIN, WORK[JVEC]);
      ENRMAX = max(ENRMAX, WORK[JVEC]);
    }
  }

  // Norm of A:

  ANORM = max(dlange(NORMA, N, N, A, LDA, WORK), UNFL);

  // Norm of E:

  ENORM = max(dlange(NORME, N, N, E, LDE, WORK), ULP);

  // Norm of error:

  // Error =  AE - EW

  dlaset('Full', N, N, ZERO, ZERO, WORK.asMatrix(N), N);

  IPAIR = 0;
  IEROW = 1;
  IECOL = 1;

  for (JCOL = 1; JCOL <= N; JCOL++) {
    if (ITRNSE == 1) {
      IEROW = JCOL;
    } else {
      IECOL = JCOL;
    }

    if (IPAIR == 0 && WI[JCOL] != ZERO) IPAIR = 1;

    if (IPAIR == 1) {
      WMAT[1][1] = WR[JCOL];
      WMAT[2][1] = -WI[JCOL];
      WMAT[1][2] = WI[JCOL];
      WMAT[2][2] = WR[JCOL];
      dgemm(TRANSE, TRANSW, N, 2, 2, ONE, E(IEROW, IECOL), LDE, WMAT, 2, ZERO,
          WORK(N * (JCOL - 1) + 1).asMatrix(N), N);
      IPAIR = 2;
    } else if (IPAIR == 2) {
      IPAIR = 0;
    } else {
      daxpy(N, WR[JCOL], E(IEROW, IECOL).asArray(), INCE,
          WORK(N * (JCOL - 1) + 1), 1);
      IPAIR = 0;
    }
  }

  dgemm(
      TRANSA, TRANSE, N, N, N, ONE, A, LDA, E, LDE, -ONE, WORK.asMatrix(N), N);

  ERRNRM = dlange('One', N, N, WORK.asMatrix(N), N, WORK(N * N + 1)) / ENORM;

  // Compute RESULT[1] (avoiding under/overflow)

  if (ANORM > ERRNRM) {
    RESULT[1] = (ERRNRM / ANORM) / ULP;
  } else {
    if (ANORM < ONE) {
      RESULT[1] = ONE / ULP;
    } else {
      RESULT[1] = min(ERRNRM / ANORM, ONE) / ULP;
    }
  }

  // Compute RESULT[2] : the normalization error in E.

  RESULT[2] = max((ENRMAX - ONE).abs(), (ENRMIN - ONE).abs()) / (N * ULP);
}
