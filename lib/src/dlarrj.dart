import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void dlarrj(
  final int N,
  final Array<double> D_,
  final Array<double> E2_,
  final int IFIRST,
  final int ILAST,
  final double RTOL,
  final int OFFSET,
  final Array<double> W_,
  final Array<double> WERR_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final double PIVMIN,
  final double SPDIAM,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E2 = E2_.having();
  final W = W_.having();
  final WERR = WERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5;
  int MAXITR;
  int CNT, I, I1, I2, II, ITER, J, K, NEXT, NINT, OLNINT, P, PREV, SAVI1;
  double DPLUS, FAC = 0, LEFT, MID, RIGHT, S, TMP, WIDTH;

  INFO.value = 0;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  MAXITR = ((log(SPDIAM + PIVMIN) - log(PIVMIN)) ~/ log(TWO)) + 2;

  // Initialize unconverged intervals in [ WORK[2*I-1], WORK[2*I] ].
  // The Sturm Count, Count( WORK[2*I-1] ) is arranged to be I-1, while
  // Count( WORK[2*I] ) is stored in IWORK[ 2*I ]. The integer IWORK[ 2*I-1 ]
  // for an unconverged interval is set to the index of the next unconverged
  // interval, and is -1 or 0 for a converged interval. Thus a linked
  // list of unconverged intervals is set up.

  I1 = IFIRST;
  I2 = ILAST;
  // The number of unconverged intervals
  NINT = 0;
  // The last unconverged interval found
  PREV = 0;
  for (I = I1; I <= I2; I++) {
    K = 2 * I;
    II = I - OFFSET;
    LEFT = W[II] - WERR[II];
    MID = W[II];
    RIGHT = W[II] + WERR[II];
    WIDTH = RIGHT - MID;
    TMP = max((LEFT).abs(), (RIGHT).abs());

    // The following test prevents the test of converged intervals
    if (WIDTH < RTOL * TMP) {
      // This interval has already converged and does not need refinement.
      // (Note that the gaps might change through refining the
      // eigenvalues, however, they can only get bigger.)
      // Remove it from the list.
      IWORK[K - 1] = -1;
      // Make sure that I1 always points to the first unconverged interval
      if ((I == I1) && (I < I2)) I1 = I + 1;
      if ((PREV >= I1) && (I <= I2)) IWORK[2 * PREV - 1] = I + 1;
    } else {
      // unconverged interval found
      PREV = I;
      // Make sure that [LEFT,RIGHT] contains the desired eigenvalue

      // Do while( CNT(LEFT) > I-1 )

      FAC = ONE;
      while (true) {
        CNT = 0;
        S = LEFT;
        DPLUS = D[1] - S;
        if (DPLUS < ZERO) CNT = CNT + 1;
        for (J = 2; J <= N; J++) {
          DPLUS = D[J] - S - E2[J - 1] / DPLUS;
          if (DPLUS < ZERO) CNT = CNT + 1;
        }
        if (CNT > I - 1) {
          LEFT = LEFT - WERR[II] * FAC;
          FAC = TWO * FAC;
          continue;
        }
        break;
      }

      // Do while( CNT(RIGHT) < I )

      FAC = ONE;
      while (true) {
        CNT = 0;
        S = RIGHT;
        DPLUS = D[1] - S;
        if (DPLUS < ZERO) CNT = CNT + 1;
        for (J = 2; J <= N; J++) {
          DPLUS = D[J] - S - E2[J - 1] / DPLUS;
          if (DPLUS < ZERO) CNT = CNT + 1;
        }
        if (CNT < I) {
          RIGHT = RIGHT + WERR[II] * FAC;
          FAC = TWO * FAC;
          continue;
        }
        break;
      }
      NINT++;
      IWORK[K - 1] = I + 1;
      IWORK[K] = CNT;
    }
    WORK[K - 1] = LEFT;
    WORK[K] = RIGHT;
  }

  SAVI1 = I1;

  // Do while( NINT > 0 ), i.e. there are still unconverged intervals
  // and while (ITER < MAXITR)

  ITER = 0;
  do {
    PREV = I1 - 1;
    I = I1;
    OLNINT = NINT;

    for (P = 1; P <= OLNINT; P++) {
      K = 2 * I;
      II = I - OFFSET;
      NEXT = IWORK[K - 1];
      LEFT = WORK[K - 1];
      RIGHT = WORK[K];
      MID = HALF * (LEFT + RIGHT);

      // semiwidth of interval
      WIDTH = RIGHT - MID;
      TMP = max((LEFT).abs(), (RIGHT).abs());
      if ((WIDTH < RTOL * TMP) || (ITER == MAXITR)) {
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

      CNT = 0;
      S = MID;
      DPLUS = D[1] - S;
      if (DPLUS < ZERO) CNT = CNT + 1;
      for (J = 2; J <= N; J++) {
        DPLUS = D[J] - S - E2[J - 1] / DPLUS;
        if (DPLUS < ZERO) CNT = CNT + 1;
      }
      if (CNT <= I - 1) {
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
  for (I = SAVI1; I <= ILAST; I++) {
    K = 2 * I;
    II = I - OFFSET;
    // All intervals marked by '0' have been refined.
    if (IWORK[K - 1] == 0) {
      W[II] = HALF * (WORK[K - 1] + WORK[K]);
      WERR[II] = WORK[K] - W[II];
    }
  }
}
