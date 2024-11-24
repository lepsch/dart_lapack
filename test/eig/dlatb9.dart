// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/lsamen.dart';

bool _FIRST = true;
double _BADC1 = 0, _BADC2 = 0, _EPS = 0, _LARGE = 0, _SMALL = 0;

({
  String TYPE,
  int KLA,
  int KUA,
  int KLB,
  int KUB,
  double ANORM,
  double BNORM,
  int MODEA,
  int MODEB,
  double CNDNMA,
  double CNDNMB,
  String DISTA,
  String DISTB,
}) dlatb9(
  final String PATH,
  final int IMAT,
  final int M,
  final int P,
  final int N,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const SHRINK = 0.25, TENTH = 0.1;
  const ONE = 1.0, TEN = 1.0e+1;

  final String TYPE, DISTA, DISTB;
  int KLA = 0, KUA = 0, KLB = 0, KUB = 0, MODEA, MODEB;
  double ANORM, BNORM, CNDNMA, CNDNMB;

  // Set some constants for use in the subroutine.

  if (_FIRST) {
    _FIRST = false;
    _EPS = dlamch('Precision');
    _BADC2 = TENTH / _EPS;
    _BADC1 = sqrt(_BADC2);
    _SMALL = dlamch('Safe minimum');
    _LARGE = ONE / _SMALL;
    _SMALL = SHRINK * (_SMALL / _EPS);
    _LARGE = ONE / _SMALL;
  }

  // Set some parameters we don't plan to change.

  TYPE = 'N';
  DISTA = 'S';
  DISTB = 'S';
  MODEA = 3;
  MODEB = 4;

  // Set the lower and upper bandwidths.

  if (lsamen(3, PATH, 'GRQ') ||
      lsamen(3, PATH, 'LSE') ||
      lsamen(3, PATH, 'GSV')) {
    // A: M by N, B: P by N

    if (IMAT == 1) {
      // A: diagonal, B: upper triangular

      KLA = 0;
      KUA = 0;
      KLB = 0;
      KUB = max(N - 1, 0);
    } else if (IMAT == 2) {
      // A: upper triangular, B: upper triangular

      KLA = 0;
      KUA = max(N - 1, 0);
      KLB = 0;
      KUB = max(N - 1, 0);
    } else if (IMAT == 3) {
      // A: lower triangular, B: upper triangular

      KLA = max(M - 1, 0);
      KUA = 0;
      KLB = 0;
      KUB = max(N - 1, 0);
    } else {
      // A: general dense, B: general dense

      KLA = max(M - 1, 0);
      KUA = max(N - 1, 0);
      KLB = max(P - 1, 0);
      KUB = max(N - 1, 0);
    }
  } else if (lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GLM')) {
    // A: N by M, B: N by P

    if (IMAT == 1) {
      // A: diagonal, B: lower triangular

      KLA = 0;
      KUA = 0;
      KLB = max(N - 1, 0);
      KUB = 0;
    } else if (IMAT == 2) {
      // A: lower triangular, B: diagonal

      KLA = max(N - 1, 0);
      KUA = 0;
      KLB = 0;
      KUB = 0;
    } else if (IMAT == 3) {
      // A: lower triangular, B: upper triangular

      KLA = max(N - 1, 0);
      KUA = 0;
      KLB = 0;
      KUB = max(P - 1, 0);
    } else {
      // A: general dense, B: general dense

      KLA = max(N - 1, 0);
      KUA = max(M - 1, 0);
      KLB = max(N - 1, 0);
      KUB = max(P - 1, 0);
    }
  }

  // Set the condition number and norm.

  CNDNMA = TEN * TEN;
  CNDNMB = TEN;
  if (lsamen(3, PATH, 'GQR') ||
      lsamen(3, PATH, 'GRQ') ||
      lsamen(3, PATH, 'GSV')) {
    if (IMAT == 5) {
      CNDNMA = _BADC1;
      CNDNMB = _BADC1;
    } else if (IMAT == 6) {
      CNDNMA = _BADC2;
      CNDNMB = _BADC2;
    } else if (IMAT == 7) {
      CNDNMA = _BADC1;
      CNDNMB = _BADC2;
    } else if (IMAT == 8) {
      CNDNMA = _BADC2;
      CNDNMB = _BADC1;
    }
  }

  ANORM = TEN;
  BNORM = TEN * TEN * TEN;
  if (lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GRQ')) {
    if (IMAT == 7) {
      ANORM = _SMALL;
      BNORM = _LARGE;
    } else if (IMAT == 8) {
      ANORM = _LARGE;
      BNORM = _SMALL;
    }
  }

  if (N <= 1) {
    CNDNMA = ONE;
    CNDNMB = ONE;
  }

  return (
    TYPE: TYPE,
    KLA: KLA,
    KUA: KUA,
    KLB: KLB,
    KUB: KUB,
    ANORM: ANORM,
    BNORM: BNORM,
    MODEA: MODEA,
    MODEB: MODEB,
    CNDNMA: CNDNMA,
    CNDNMB: CNDNMB,
    DISTA: DISTA,
    DISTB: DISTB,
  );
}
