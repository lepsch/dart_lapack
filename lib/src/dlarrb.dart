// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlaneg.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarrb(
  final int N,
  final Array<double> D_,
  final Array<double> LLD_,
  final int IFIRST,
  final int ILAST,
  final double RTOL1,
  final double RTOL2,
  final int OFFSET,
  final Array<double> W_,
  final Array<double> WGAP_,
  final Array<double> WERR_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final double PIVMIN,
  final double SPDIAM,
  final int TWIST,
  final Box<int> INFO,
) {
  final D = D_.having();
  final LLD = LLD_.having();
  final W = W_.having();
  final WGAP = WGAP_.having();
  final WERR = WERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, TWO = 2.0, HALF = 0.5;
  int MAXITR;
  int I, I1, II, IP, ITER, K, NEGCNT, NEXT, NINT, OLNINT, PREV, R;
  double BACK = 0,
      CVRGD,
      GAP,
      LEFT = 0,
      LGAP,
      MID,
      MNWDTH,
      RGAP,
      RIGHT,
      TMP,
      WIDTH;

  INFO.value = 0;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  MAXITR = (log(SPDIAM + PIVMIN) - log(PIVMIN)) ~/ log(TWO) + 2;
  MNWDTH = TWO * PIVMIN;

  R = TWIST;
  if ((R < 1) || (R > N)) R = N;

  // Initialize unconverged intervals in [ WORK[2*I-1], WORK[2*I] ].
  // The Sturm Count, Count( WORK[2*I-1] ) is arranged to be I-1, while
  // Count( WORK[2*I] ) is stored in IWORK[ 2*I ]. The integer IWORK[ 2*I-1 ]
  // for an unconverged interval is set to the index of the next unconverged
  // interval, and is -1 or 0 for a converged interval. Thus a linked
  // list of unconverged intervals is set up.

  I1 = IFIRST;
  // The number of unconverged intervals
  NINT = 0;
  // The last unconverged interval found
  PREV = 0;

  RGAP = WGAP[I1 - OFFSET];
  for (I = I1; I <= ILAST; I++) {
    K = 2 * I;
    II = I - OFFSET;
    LEFT = W[II] - WERR[II];
    RIGHT = W[II] + WERR[II];
    LGAP = RGAP;
    RGAP = WGAP[II];
    GAP = min(LGAP, RGAP);

    // Make sure that [LEFT,RIGHT] contains the desired eigenvalue
    // Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT

    // Do while( NEGCNT(LEFT) > I-1 )

    BACK = WERR[II];
    while (true) {
      NEGCNT = dlaneg(N, D, LLD, LEFT, PIVMIN, R);
      if (NEGCNT <= I - 1) break;
      LEFT -= BACK;
      BACK = TWO * BACK;
    }

    // Do while( NEGCNT(RIGHT) < I )
    // Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT

    BACK = WERR[II];
    while (true) {
      NEGCNT = dlaneg(N, D, LLD, RIGHT, PIVMIN, R);
      if (NEGCNT >= I) break;
      RIGHT += BACK;
      BACK = TWO * BACK;
    }
    WIDTH = HALF * (LEFT - RIGHT).abs();
    TMP = max(LEFT.abs(), RIGHT.abs());
    CVRGD = max(RTOL1 * GAP, RTOL2 * TMP);
    if (WIDTH <= CVRGD || WIDTH <= MNWDTH) {
      // This interval has already converged and does not need refinement.
      // (Note that the gaps might change through refining the
      // eigenvalues, however, they can only get bigger.)
      // Remove it from the list.
      IWORK[K - 1] = -1;
      // Make sure that I1 always points to the first unconverged interval
      if ((I == I1) && (I < ILAST)) I1 = I + 1;
      if ((PREV >= I1) && (I <= ILAST)) IWORK[2 * PREV - 1] = I + 1;
    } else {
      // unconverged interval found
      PREV = I;
      NINT++;
      IWORK[K - 1] = I + 1;
      IWORK[K] = NEGCNT;
    }
    WORK[K - 1] = LEFT;
    WORK[K] = RIGHT;
  }

  // Do while( NINT > 0 ), i.e. there are still unconverged intervals
  // and while (ITER < MAXITR)

  ITER = 0;
  do {
    PREV = I1 - 1;
    I = I1;
    OLNINT = NINT;

    for (IP = 1; IP <= OLNINT; IP++) {
      K = 2 * I;
      II = I - OFFSET;
      RGAP = WGAP[II];
      LGAP = RGAP;
      if (II > 1) LGAP = WGAP[II - 1];
      GAP = min(LGAP, RGAP);
      NEXT = IWORK[K - 1];
      LEFT = WORK[K - 1];
      RIGHT = WORK[K];
      MID = HALF * (LEFT + RIGHT);

      // semiwidth of interval
      WIDTH = RIGHT - MID;
      TMP = max(LEFT.abs(), RIGHT.abs());
      CVRGD = max(RTOL1 * GAP, RTOL2 * TMP);
      if ((WIDTH <= CVRGD) || (WIDTH <= MNWDTH) || (ITER == MAXITR)) {
        // reduce number of unconverged intervals
        NINT--;
        // Mark interval as converged.
        IWORK[K - 1] = 0;
        if (I1 == I) {
          I1 = NEXT;
        } else {
          // Prev holds the last unconverged interval previously examined
          if (PREV >= I1) IWORK[2 * PREV - 1] = NEXT;
        }
        I = NEXT;
        continue;
      }
      PREV = I;

      // Perform one bisection step

      NEGCNT = dlaneg(N, D, LLD, MID, PIVMIN, R);
      if (NEGCNT <= I - 1) {
        WORK[K - 1] = MID;
      } else {
        WORK[K] = MID;
      }
      I = NEXT;
    }
    ITER++;
    // do another loop if there are still unconverged intervals
    // However, in the last iteration, all intervals are accepted
    // since this is the best we can do.
  } while ((NINT > 0) && (ITER <= MAXITR));

  // At this point, all the intervals have converged
  for (I = IFIRST; I <= ILAST; I++) {
    K = 2 * I;
    II = I - OFFSET;
    // All intervals marked by '0' have been refined.
    if (IWORK[K - 1] == 0) {
      W[II] = HALF * (WORK[K - 1] + WORK[K]);
      WERR[II] = WORK[K] - W[II];
    }
  }

  for (I = IFIRST + 1; I <= ILAST; I++) {
    K = 2 * I;
    II = I - OFFSET;
    WGAP[II - 1] = max(ZERO, W[II] - WERR[II] - W[II - 1] - WERR[II - 1]);
  }
}
