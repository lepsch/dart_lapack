// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/matrix.dart';

int dlaneg(
  final int N,
  final Array<double> D_,
  final Array<double> LLD_,
  final double SIGMA,
  final double PIVMIN,
  final int R,
) {
  final D = D_.having();
  final LLD = LLD_.having();
  const ZERO = 0.0, ONE = 1.0;
  // Some architectures propagate Infinities and NaNs very slowly, so
  // the code computes counts in BLKLEN chunks.  Then a NaN can
  // propagate at most BLKLEN columns before being detected.  This is
  // not a general tuning parameter; it needs only to be just large
  // enough that the overhead is tiny in common cases.
  const BLKLEN = 128;
  int BJ, J, NEG1, NEG2, NEGCNT;
  double BSAV, DMINUS, DPLUS, GAMMA, P, T, TMP;
  bool SAWNAN;

  NEGCNT = 0;

  // I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
  T = -SIGMA;
  for (BJ = 1; BJ <= R - 1; BJ += BLKLEN) {
    NEG1 = 0;
    BSAV = T;
    for (J = BJ; J <= min(BJ + BLKLEN - 1, R - 1); J++) {
      DPLUS = D[J] + T;
      if (DPLUS < ZERO) NEG1++;
      TMP = T / DPLUS;
      T = TMP * LLD[J] - SIGMA;
    }
    SAWNAN = disnan(T);
    // Run a slower version of the above loop if a NaN is detected.
    // A NaN should occur only with a zero pivot after an infinite
    // pivot.  In that case, substituting 1 for T/DPLUS is the
    // correct limit.
    if (SAWNAN) {
      NEG1 = 0;
      T = BSAV;
      for (J = BJ; J <= min(BJ + BLKLEN - 1, R - 1); J++) {
        DPLUS = D[J] + T;
        if (DPLUS < ZERO) NEG1++;
        TMP = T / DPLUS;
        if (disnan(TMP)) TMP = ONE;
        T = TMP * LLD[J] - SIGMA;
      }
    }
    NEGCNT += NEG1;
  }

  // II) lower part: L D L^T - SIGMA I = U- D- U-^T
  P = D[N] - SIGMA;
  for (BJ = N - 1; BJ >= R; BJ -= BLKLEN) {
    NEG2 = 0;
    BSAV = P;
    for (J = BJ; J >= max(BJ - BLKLEN + 1, R); J--) {
      DMINUS = LLD[J] + P;
      if (DMINUS < ZERO) NEG2++;
      TMP = P / DMINUS;
      P = TMP * D[J] - SIGMA;
    }
    SAWNAN = disnan(P);
    // As above, run a slower version that substitutes 1 for Inf/Inf.

    if (SAWNAN) {
      NEG2 = 0;
      P = BSAV;
      for (J = BJ; J >= max(BJ - BLKLEN + 1, R); J--) {
        DMINUS = LLD[J] + P;
        if (DMINUS < ZERO) NEG2++;
        TMP = P / DMINUS;
        if (disnan(TMP)) TMP = ONE;
        P = TMP * D[J] - SIGMA;
      }
    }
    NEGCNT += NEG2;
  }

  // III) Twist index
  //   T was shifted by SIGMA initially.
  GAMMA = (T + SIGMA) + P;
  if (GAMMA < ZERO) NEGCNT++;

  return NEGCNT;
}
