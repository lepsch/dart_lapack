// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dget52(
  final bool LEFT,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> E_,
  final int LDE,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Array<double> WORK_,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final E = E_.having(ld: LDE);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  bool ILCPLX;
  String NORMAB, TRANS;
  int J, JVEC;
  double ABMAX,
      ACOEF,
      ALFMAX,
      ANORM,
      BCOEFI,
      BCOEFR,
      BETMAX,
      BNORM,
      ENORM,
      ENRMER,
      ERRNRM,
      SAFMAX,
      SAFMIN,
      SALFI,
      SALFR,
      SBETA,
      SCALE,
      TEMP1,
      ULP;

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  SAFMIN = dlamch('Safe minimum');
  SAFMAX = ONE / SAFMIN;
  ULP = dlamch('Epsilon') * dlamch('Base');

  if (LEFT) {
    TRANS = 'T';
    NORMAB = 'I';
  } else {
    TRANS = 'N';
    NORMAB = 'O';
  }

  // Norm of A, B, and E:

  ANORM = max(dlange(NORMAB, N, N, A, LDA, WORK), SAFMIN);
  BNORM = max(dlange(NORMAB, N, N, B, LDB, WORK), SAFMIN);
  ENORM = max(dlange('O', N, N, E, LDE, WORK), ULP);
  ALFMAX = SAFMAX / max(ONE, BNORM);
  BETMAX = SAFMAX / max(ONE, ANORM);

  // Compute error matrix.
  // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

  ILCPLX = false;
  for (JVEC = 1; JVEC <= N; JVEC++) {
    if (ILCPLX) {
      // 2nd Eigenvalue/-vector of pair -- do nothing

      ILCPLX = false;
    } else {
      SALFR = ALPHAR[JVEC];
      SALFI = ALPHAI[JVEC];
      SBETA = BETA[JVEC];
      if (SALFI == ZERO) {
        // Real eigenvalue and -vector

        ABMAX = max(SALFR.abs(), SBETA.abs());
        if (SALFR.abs() > ALFMAX || SBETA.abs() > BETMAX || ABMAX < ONE) {
          SCALE = ONE / max(ABMAX, SAFMIN);
          SALFR = SCALE * SALFR;
          SBETA = SCALE * SBETA;
        }
        SCALE =
            ONE / max(SALFR.abs() * BNORM, max(SBETA.abs() * ANORM, SAFMIN));
        ACOEF = SCALE * SBETA;
        BCOEFR = SCALE * SALFR;
        dgemv(TRANS, N, N, ACOEF, A, LDA, E(1, JVEC).asArray(), 1, ZERO,
            WORK(N * (JVEC - 1) + 1), 1);
        dgemv(TRANS, N, N, -BCOEFR, B, LDA, E(1, JVEC).asArray(), 1, ONE,
            WORK(N * (JVEC - 1) + 1), 1);
      } else {
        // Complex conjugate pair

        ILCPLX = true;
        if (JVEC == N) {
          RESULT[1] = TEN / ULP;
          return;
        }
        ABMAX = max(SALFR.abs() + SALFI.abs(), SBETA.abs());
        if (SALFR.abs() + SALFI.abs() > ALFMAX ||
            SBETA.abs() > BETMAX ||
            ABMAX < ONE) {
          SCALE = ONE / max(ABMAX, SAFMIN);
          SALFR = SCALE * SALFR;
          SALFI = SCALE * SALFI;
          SBETA = SCALE * SBETA;
        }
        SCALE = ONE /
            max((SALFR.abs() + SALFI.abs()) * BNORM,
                max(SBETA.abs() * ANORM, SAFMIN));
        ACOEF = SCALE * SBETA;
        BCOEFR = SCALE * SALFR;
        BCOEFI = SCALE * SALFI;
        if (LEFT) {
          BCOEFI = -BCOEFI;
        }

        dgemv(TRANS, N, N, ACOEF, A, LDA, E(1, JVEC).asArray(), 1, ZERO,
            WORK(N * (JVEC - 1) + 1), 1);
        dgemv(TRANS, N, N, -BCOEFR, B, LDA, E(1, JVEC).asArray(), 1, ONE,
            WORK(N * (JVEC - 1) + 1), 1);
        dgemv(TRANS, N, N, BCOEFI, B, LDA, E(1, JVEC + 1).asArray(), 1, ONE,
            WORK(N * (JVEC - 1) + 1), 1);

        dgemv(TRANS, N, N, ACOEF, A, LDA, E(1, JVEC + 1).asArray(), 1, ZERO,
            WORK(N * JVEC + 1), 1);
        dgemv(TRANS, N, N, -BCOEFI, B, LDA, E(1, JVEC).asArray(), 1, ONE,
            WORK(N * JVEC + 1), 1);
        dgemv(TRANS, N, N, -BCOEFR, B, LDA, E(1, JVEC + 1).asArray(), 1, ONE,
            WORK(N * JVEC + 1), 1);
      }
    }
  }

  ERRNRM =
      dlange('One', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1)) /
          ENORM;

  // Compute RESULT[1]

  RESULT[1] = ERRNRM / ULP;

  // Normalization of E:

  ENRMER = ZERO;
  ILCPLX = false;
  for (JVEC = 1; JVEC <= N; JVEC++) {
    if (ILCPLX) {
      ILCPLX = false;
    } else {
      TEMP1 = ZERO;
      if (ALPHAI[JVEC] == ZERO) {
        for (J = 1; J <= N; J++) {
          TEMP1 = max(TEMP1, E[J][JVEC].abs());
        }
        ENRMER = max(ENRMER, (TEMP1 - ONE).abs());
      } else {
        ILCPLX = true;
        for (J = 1; J <= N; J++) {
          TEMP1 = max(TEMP1, E[J][JVEC].abs() + E[J][JVEC + 1].abs());
        }
        ENRMER = max(ENRMER, (TEMP1 - ONE).abs());
      }
    }
  }

  // Compute RESULT[2] : the normalization error in E.

  RESULT[2] = ENRMER / (N * ULP);
}
