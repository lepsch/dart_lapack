import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlarrf(
  final int N,
  final Array<double> D_,
  final Array<double> L_,
  final Array<double> LD_,
  final int CLSTRT,
  final int CLEND,
  final Array<double> W_,
  final Array<double> WGAP_,
  final Array<double> WERR_,
  final double SPDIAM,
  final double CLGAPL,
  final double CLGAPR,
  final double PIVMIN,
  final Box<double> SIGMA,
  final Array<double> DPLUS_,
  final Array<double> LPLUS_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final L = L_.having();
  final LD = LD_.having();
  final W = W_.having();
  final WGAP = WGAP_.having();
  final WERR = WERR_.having();
  final DPLUS = DPLUS_.having();
  final LPLUS = LPLUS_.having();
  final WORK = WORK_.having();
  const ONE = 1.0,
      TWO = 2.0,
      FOUR = 4.0,
      QUART = 0.25,
      MAXGROWTH1 = 8.0,
      MAXGROWTH2 = 8.0;
  bool DORRR1, FORCER, NOFAIL, SAWNAN1, SAWNAN2, TRYRRR1;
  int I, INDX = 0, KTRY, SHIFT;
  const KTRYMAX = 1, SLEFT = 1, SRIGHT = 2;
  double AVGAP,
      BESTSHIFT,
      CLWDTH,
      EPS,
      FACT,
      FAIL,
      FAIL2,
      GROWTHBOUND,
      LDELTA,
      LDMAX,
      LSIGMA,
      MAX1,
      MAX2,
      MINGAP,
      OLDP,
      PROD,
      RDELTA,
      RDMAX,
      RRR1,
      RRR2,
      RSIGMA,
      S,
      SMLGROWTH,
      TMP,
      ZNM2;

  INFO.value = 0;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  FACT = pow(2, KTRYMAX).toDouble();
  EPS = dlamch('Precision');
  SHIFT = 0;
  FORCER = false;

  // Note that we cannot guarantee that for any of the shifts tried,
  // the factorization has a small or even moderate element growth.
  // There could be Ritz values at both ends of the cluster and despite
  // backing off, there are examples where all factorizations tried
  // (in IEEE mode, allowing zero pivots & infinities) have INFINITE
  // element growth.
  // For this reason, we should use PIVMIN in this subroutine so that at
  // least the L D L^T factorization exists. It can be checked afterwards
  // whether the element growth caused bad residuals/orthogonality.

  // Decide whether the code should accept the best among all
  // representations despite large element growth or signal INFO.value=1
  // Setting NOFAIL to false for quick fix for bug 113
  NOFAIL = false;

  // Compute the average gap length of the cluster
  CLWDTH = (W[CLEND] - W[CLSTRT]).abs() + WERR[CLEND] + WERR[CLSTRT];
  AVGAP = CLWDTH / (CLEND - CLSTRT).toDouble();
  MINGAP = min(CLGAPL, CLGAPR);
  // Initial values for shifts to both ends of cluster
  LSIGMA = min(W[CLSTRT], W[CLEND]) - WERR[CLSTRT];
  RSIGMA = max(W[CLSTRT], W[CLEND]) + WERR[CLEND];

  // Use a small fudge to make sure that we really shift to the outside
  LSIGMA -= (LSIGMA).abs() * FOUR * EPS;
  RSIGMA += (RSIGMA).abs() * FOUR * EPS;

  // Compute upper bounds for how much to back off the initial shifts
  LDMAX = QUART * MINGAP + TWO * PIVMIN;
  RDMAX = QUART * MINGAP + TWO * PIVMIN;

  LDELTA = max(AVGAP, WGAP[CLSTRT]) / FACT;
  RDELTA = max(AVGAP, WGAP[CLEND - 1]) / FACT;

  // Initialize the record of the best representation found

  S = dlamch('S');
  SMLGROWTH = ONE / S;
  FAIL = (N - 1).toDouble() * MINGAP / (SPDIAM * EPS);
  FAIL2 = (N - 1).toDouble() * MINGAP / (SPDIAM * sqrt(EPS));
  BESTSHIFT = LSIGMA;

  // while (KTRY <= KTRYMAX)
  KTRY = 0;
  GROWTHBOUND = MAXGROWTH1 * SPDIAM;

  var ineffective = false;
  do {
    SAWNAN1 = false;
    SAWNAN2 = false;
    // Ensure that we do not back off too much of the initial shifts
    LDELTA = min(LDMAX, LDELTA);
    RDELTA = min(RDMAX, RDELTA);

    // Compute the element growth when shifting to both ends of the cluster
    // accept the shift if there is no element growth at one of the two ends

    // Left end
    S = -LSIGMA;
    DPLUS[1] = D[1] + S;
    if ((DPLUS[1]).abs() < PIVMIN) {
      DPLUS[1] = -PIVMIN;
      // Need to set SAWNAN1 because refined RRR test should not be used
      // in this case
      SAWNAN1 = true;
    }
    MAX1 = (DPLUS[1]).abs();
    for (I = 1; I <= N - 1; I++) {
      LPLUS[I] = LD[I] / DPLUS[I];
      S = S * LPLUS[I] * L[I] - LSIGMA;
      DPLUS[I + 1] = D[I + 1] + S;
      if ((DPLUS[I + 1]).abs() < PIVMIN) {
        DPLUS[I + 1] = -PIVMIN;
        // Need to set SAWNAN1 because refined RRR test should not be used
        // in this case
        SAWNAN1 = true;
      }
      MAX1 = max(MAX1, (DPLUS[I + 1]).abs());
    }
    SAWNAN1 = SAWNAN1 || disnan(MAX1);
    if (FORCER || (MAX1 <= GROWTHBOUND && !SAWNAN1)) {
      SIGMA.value = LSIGMA;
      SHIFT = SLEFT;
      break;
    }

    // Right end
    S = -RSIGMA;
    WORK[1] = D[1] + S;
    if ((WORK[1]).abs() < PIVMIN) {
      WORK[1] = -PIVMIN;
      // Need to set SAWNAN2 because refined RRR test should not be used
      // in this case
      SAWNAN2 = true;
    }
    MAX2 = (WORK[1]).abs();
    for (I = 1; I <= N - 1; I++) {
      WORK[N + I] = LD[I] / WORK[I];
      S = S * WORK[N + I] * L[I] - RSIGMA;
      WORK[I + 1] = D[I + 1] + S;
      if ((WORK[I + 1]).abs() < PIVMIN) {
        WORK[I + 1] = -PIVMIN;
        // Need to set SAWNAN2 because refined RRR test should not be used
        // in this case
        SAWNAN2 = true;
      }
      MAX2 = max(MAX2, (WORK[I + 1]).abs());
    }
    SAWNAN2 = SAWNAN2 || disnan(MAX2);
    if (FORCER || (MAX2 <= GROWTHBOUND && !SAWNAN2)) {
      SIGMA.value = RSIGMA;
      SHIFT = SRIGHT;
      break;
    }
    // If we are at this point, both shifts led to too much element growth

    // Record the better of the two shifts (provided it didn't lead to NaN)
    // both MAX1 and MAX2 are NaN if SAWNAN1 && SAWNAN2
    if (!SAWNAN1 || !SAWNAN2) {
      if (!SAWNAN1) {
        INDX = 1;
        if (MAX1 <= SMLGROWTH) {
          SMLGROWTH = MAX1;
          BESTSHIFT = LSIGMA;
        }
      }
      if (!SAWNAN2) {
        if (SAWNAN1 || MAX2 <= MAX1) INDX = 2;
        if (MAX2 <= SMLGROWTH) {
          SMLGROWTH = MAX2;
          BESTSHIFT = RSIGMA;
        }
      }

      // If we are here, both the left and the right shift led to
      // element growth. If the element growth is moderate, then
      // we may still accept the representation, if it passes a
      // refined test for RRR. This test supposes that no NaN occurred.
      // Moreover, we use the refined RRR test only for isolated clusters.
      if ((CLWDTH < MINGAP / 128.toDouble()) &&
          (min(MAX1, MAX2) < FAIL2) &&
          (!SAWNAN1) &&
          (!SAWNAN2)) {
        DORRR1 = true;
      } else {
        DORRR1 = false;
      }
      TRYRRR1 = true;
      if (TRYRRR1 && DORRR1) {
        if (INDX == 1) {
          TMP = (DPLUS[N]).abs();
          ZNM2 = ONE;
          PROD = ONE;
          OLDP = ONE;
          for (I = N - 1; I >= 1; I--) {
            if (PROD <= EPS) {
              PROD = ((DPLUS[I + 1] * WORK[N + I + 1]) /
                      (DPLUS[I] * WORK[N + I])) *
                  OLDP;
            } else {
              PROD = PROD * (WORK[N + I]).abs();
            }
            OLDP = PROD;
            ZNM2 += pow(PROD, 2);
            TMP = max(TMP, (DPLUS[I] * PROD).abs());
          }
          RRR1 = TMP / (SPDIAM * sqrt(ZNM2));
          if (RRR1 <= MAXGROWTH2) {
            SIGMA.value = LSIGMA;
            SHIFT = SLEFT;
            break;
          }
        } else if (INDX == 2) {
          TMP = (WORK[N]).abs();
          ZNM2 = ONE;
          PROD = ONE;
          OLDP = ONE;
          for (I = N - 1; I >= 1; I--) {
            if (PROD <= EPS) {
              PROD =
                  ((WORK[I + 1] * LPLUS[I + 1]) / (WORK[I] * LPLUS[I])) * OLDP;
            } else {
              PROD = PROD * (LPLUS[I]).abs();
            }
            OLDP = PROD;
            ZNM2 += pow(PROD, 2);
            TMP = max(TMP, (WORK[I] * PROD).abs());
          }
          RRR2 = TMP / (SPDIAM * sqrt(ZNM2));
          if (RRR2 <= MAXGROWTH2) {
            SIGMA.value = RSIGMA;
            SHIFT = SRIGHT;
            break;
          }
        }
      }
    }

    ineffective = KTRY < KTRYMAX;
    if (ineffective) {
      // If we are here, both shifts ineffective also the RRR test.
      // Back off to the outside
      LSIGMA = max(LSIGMA - LDELTA, LSIGMA - LDMAX);
      RSIGMA = min(RSIGMA + RDELTA, RSIGMA + RDMAX);
      LDELTA = TWO * LDELTA;
      RDELTA = TWO * RDELTA;
      KTRY++;
    } else {
      // None of the representations investigated satisfied our
      // criteria. Take the best one we found.
      if ((SMLGROWTH < FAIL) || NOFAIL) {
        LSIGMA = BESTSHIFT;
        RSIGMA = BESTSHIFT;
        FORCER = true;
        ineffective = true;
      } else {
        INFO.value = 1;
        return;
      }
    }
  } while (ineffective);

  if (SHIFT == SLEFT) {
  } else if (SHIFT == SRIGHT) {
    // store new L and D back into DPLUS, LPLUS
    dcopy(N, WORK, 1, DPLUS, 1);
    dcopy(N - 1, WORK(N + 1), 1, LPLUS, 1);
  }
}
