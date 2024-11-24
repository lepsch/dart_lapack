// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'dsvdct.dart';

void dsvdch(
  final int N,
  final Array<double> S_,
  final Array<double> E_,
  final Array<double> SVD_,
  final double TOL,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final S = S_.having();
  final E = E_.having();
  final SVD = SVD_.having();
  const ONE = 1.0;
  const ZERO = 0.0;
  int BPNT, COUNT, TPNT;
  double EPS, LOWER, OVFL, TUPPR, UNFL, UNFLEP, UPPER;
  final NUML = Box(0), NUMU = Box(0);

  // Get machine constants

  INFO.value = 0;
  if (N <= 0) return;
  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');
  EPS = dlamch('Epsilon') * dlamch('Base');

  // UNFLEP is chosen so that when an eigenvalue is multiplied by the
  // scale factor sqrt(OVFL)*sqrt(sqrt(UNFL))/MX in DSVDCT, it exceeds
  // sqrt(UNFL), which is the lower limit for DSVDCT.

  UNFLEP = (sqrt(sqrt(UNFL)) / sqrt(OVFL)) * SVD[1] + UNFL / EPS;

  // The value of EPS works best when TOL >= 10.

  EPS = TOL * max(N ~/ 10, 1) * EPS;

  // TPNT points to singular value at right endpoint of interval
  // BPNT points to singular value at left  endpoint of interval

  TPNT = 1;
  BPNT = 1;

  // Begin loop over all intervals
  do {
    UPPER = (ONE + EPS) * SVD[TPNT] + UNFLEP;
    LOWER = (ONE - EPS) * SVD[BPNT] - UNFLEP;
    if (LOWER <= UNFLEP) LOWER = -UPPER;

    // Begin loop merging overlapping intervals

    while (true) {
      if (BPNT == N) break;
      TUPPR = (ONE + EPS) * SVD[BPNT + 1] + UNFLEP;
      if (TUPPR < LOWER) break;

      // Merge

      BPNT++;
      LOWER = (ONE - EPS) * SVD[BPNT] - UNFLEP;
      if (LOWER <= UNFLEP) LOWER = -UPPER;
    }

    // Count singular values in interval [ LOWER, UPPER ]

    dsvdct(N, S, E, LOWER, NUML);
    dsvdct(N, S, E, UPPER, NUMU);
    COUNT = NUMU.value - NUML.value;
    if (LOWER < ZERO) COUNT = COUNT ~/ 2;
    if (COUNT != BPNT - TPNT + 1) {
      // Wrong number of singular values in interval

      INFO.value = TPNT;
      return;
    }
    TPNT = BPNT + 1;
    BPNT = TPNT;
  } while (TPNT <= N);
}
