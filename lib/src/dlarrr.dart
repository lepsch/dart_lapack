// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarrr(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0, RELCOND = 0.999;
  int I;
  bool YESREL;
  double EPS, SAFMIN, SMLNUM, RMIN, TMP, TMP2, OFFDIG, OFFDIG2;

  // Quick return if possible

  if (N <= 0) {
    INFO.value = 0;
    return;
  }

  // As a default, do NOT go for relative-accuracy preserving computations.
  INFO.value = 1;

  SAFMIN = dlamch('Safe minimum');
  EPS = dlamch('Precision');
  SMLNUM = SAFMIN / EPS;
  RMIN = sqrt(SMLNUM);

  // Tests for relative accuracy

  // Test for scaled diagonal dominance
  // Scale the diagonal entries to one and check whether the sum of the
  // off-diagonals is less than one

  // The sdd relative error bounds have a 1/(1- 2*x) factor in them,
  // x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative
  // accuracy is promised.  In the notation of the code fragment below,
  // 1/(1 - (OFFDIG + OFFDIG2)) is the condition number.
  // We don't think it is worth going into "sdd mode" unless the relative
  // condition number is reasonable, not 1/macheps.
  // The threshold should be compatible with other thresholds used in the
  // code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds
  // to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000
  // instead of the current OFFDIG + OFFDIG2 < 1

  YESREL = true;
  OFFDIG = ZERO;
  TMP = sqrt(D[1].abs());
  if (TMP < RMIN) YESREL = false;
  if (YESREL) {
    for (I = 2; I <= N; I++) {
      TMP2 = sqrt(D[I].abs());
      if (TMP2 < RMIN) YESREL = false;
      if (!YESREL) break;
      OFFDIG2 = E[I - 1].abs() / (TMP * TMP2);
      if (OFFDIG + OFFDIG2 >= RELCOND) YESREL = false;
      if (!YESREL) break;
      TMP = TMP2;
      OFFDIG = OFFDIG2;
    }
  }

  if (YESREL) {
    INFO.value = 0;
    return;
  } else {}

  // *** MORE TO BE IMPLEMENTED ***

  // Test if the lower bidiagonal matrix L from T = L D L^T
  // (zero shift facto) is well conditioned

  // Test if the upper bidiagonal matrix U from T = U D U^T
  // (zero shift facto) is well conditioned.
  // In this case, the matrix needs to be flipped and, at the end
  // of the eigenvector computation, the flip needs to be applied
  // to the computed eigenvectors (and the support)
}
