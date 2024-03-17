import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/dlarra.dart';
import 'package:lapack/src/dlarrb.dart';
import 'package:lapack/src/dlarrc.dart';
import 'package:lapack/src/dlarrd.dart';
import 'package:lapack/src/dlarrk.dart';
import 'package:lapack/src/dlasq2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlarre(
  final String RANGE,
  final int N,
  final Box<double> VL,
  final Box<double> VU,
  final int IL,
  final int IU,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> E2_,
  final double RTOL1,
  final double RTOL2,
  final double SPLTOL,
  final Box<int> NSPLIT,
  final Array<int> ISPLIT_,
  final Box<int> M,
  final Array<double> W_,
  final Array<double> WERR_,
  final Array<double> WGAP_,
  final Array<int> IBLOCK_,
  final Array<int> INDEXW_,
  final Array<double> GERS_,
  final Box<double> PIVMIN,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final E2 = E2_.having();
  final ISPLIT = ISPLIT_.having();
  final W = W_.having();
  final WERR = WERR_.having();
  final WGAP = WGAP_.having();
  final IBLOCK = IBLOCK_.having();
  final INDEXW = INDEXW_.having();
  final GERS = GERS_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0,
      ONE = 1.0,
      TWO = 2.0,
      FOUR = 4.0,
      HNDRD = 100.0,
      PERT = 8.0,
      HALF = ONE / TWO,
      FOURTH = ONE / FOUR,
      FAC = HALF,
      MAXGROWTH = 64.0,
      FUDGE = 2.0;
  const MAXTRY = 6, ALLRNG = 1, INDRNG = 2, VALRNG = 3;
  bool FORCEB, NOREP, USEDQD;
  int I,
      IBEGIN,
      IDUM,
      IEND = 0,
      IN = 0,
      INDL = 0,
      INDU = 0,
      IRANGE = 0,
      J,
      JBLK,
      MB = 0,
      WBEGIN,
      WEND = 0;
  double AVGAP,
      BSRTOL,
      CLWDTH,
      DMAX,
      DPIVOT,
      EABS,
      EMAX,
      EOLD,
      EPS,
      GL,
      GU,
      ISLEFT,
      ISRGHT,
      RTL,
      RTOL,
      S1,
      S2,
      SAFMIN,
      SGNDEF,
      SIGMA,
      SPDIAM,
      TAU;
  final ISEED = Array<int>(4);
  final CNT = Box(0), CNT1 = Box(0), CNT2 = Box(0), IINFO = Box(0), MM = Box(0);
  final TMP = Box(0.0), TMP1 = Box(0.0);

  INFO.value = 0;
  NSPLIT.value = 0;
  M.value = 0;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  // Decode RANGE

  if (lsame(RANGE, 'A')) {
    IRANGE = ALLRNG;
  } else if (lsame(RANGE, 'V')) {
    IRANGE = VALRNG;
  } else if (lsame(RANGE, 'I')) {
    IRANGE = INDRNG;
  }

  // Get machine constants
  SAFMIN = dlamch('S');
  EPS = dlamch('P');

  // Set parameters
  RTL = sqrt(EPS);
  BSRTOL = sqrt(EPS);

  // Treat case of 1x1 matrix for quick return;
  if (N == 1) {
    if ((IRANGE == ALLRNG) ||
        ((IRANGE == VALRNG) && (D[1] > VL.value) && (D[1] <= VU.value)) ||
        ((IRANGE == INDRNG) && (IL == 1) && (IU == 1))) {
      M.value = 1;
      W[1] = D[1];
      // The computation error of the eigenvalue is zero
      WERR[1] = ZERO;
      WGAP[1] = ZERO;
      IBLOCK[1] = 1;
      INDEXW[1] = 1;
      GERS[1] = D[1];
      GERS[2] = D[1];
    }
    // store the shift for the initial RRR, which is zero in this case
    E[1] = ZERO;
    return;
  }

  // General case: tridiagonal matrix of order > 1

  // Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
  // Compute maximum off-diagonal entry and pivmin.
  GL = D[1];
  GU = D[1];
  EOLD = ZERO;
  EMAX = ZERO;
  E[N] = ZERO;
  for (I = 1; I <= N; I++) {
    WERR[I] = ZERO;
    WGAP[I] = ZERO;
    EABS = (E[I]).abs();
    if (EABS >= EMAX) {
      EMAX = EABS;
    }
    TMP1.value = EABS + EOLD;
    GERS[2 * I - 1] = D[I] - TMP1.value;
    GL = min(GL, GERS[2 * I - 1]);
    GERS[2 * I] = D[I] + TMP1.value;
    GU = max(GU, GERS[2 * I]);
    EOLD = EABS;
  }
  // The minimum pivot allowed in the Sturm sequence for T
  PIVMIN.value = SAFMIN * max(ONE, pow(EMAX, 2));
  // Compute spectral diameter. The Gerschgorin bounds give an
  // estimate that is wrong by at most a factor of sqrt(2)
  SPDIAM = GU - GL;

  // Compute splitting points
  dlarra(N, D, E, E2, SPLTOL, SPDIAM, NSPLIT, ISPLIT, IINFO);

  // Can force use of bisection instead of faster DQDS.
  // Option left in the code for future multisection work.
  FORCEB = false;

  // Initialize USEDQD, DQDS should be used for ALLRNG unless someone
  // explicitly wants bisection.
  USEDQD = ((IRANGE == ALLRNG) && (!FORCEB));

  if ((IRANGE == ALLRNG) && (!FORCEB)) {
    // Set interval [VL.value,VU.value] that contains all eigenvalues
    VL.value = GL;
    VU.value = GU;
  } else {
    // We call DLARRD to find crude approximations to the eigenvalues
    // in the desired range. In case IRANGE = INDRNG, we also obtain the
    // interval (VL.value,VU.value] that contains all the wanted eigenvalues.
    // An interval [LEFT,RIGHT] has converged if
    // RIGHT-LEFT < RTOL*max(ABS(LEFT),ABS(RIGHT))
    // DLARRD needs a WORK of size 4*N, IWORK of size 3*N
    dlarrd(
        RANGE,
        'B',
        N,
        VL.value,
        VU.value,
        IL,
        IU,
        GERS,
        BSRTOL,
        D,
        E,
        E2,
        PIVMIN.value,
        NSPLIT.value,
        ISPLIT,
        MM,
        W,
        WERR,
        VL,
        VU,
        IBLOCK,
        INDEXW,
        WORK,
        IWORK,
        IINFO);
    if (IINFO.value != 0) {
      INFO.value = -1;
      return;
    }
    // Make sure that the entries M.value+1 to N in W, WERR, IBLOCK, INDEXW are 0
    for (I = MM.value + 1; I <= N; I++) {
      W[I] = ZERO;
      WERR[I] = ZERO;
      IBLOCK[I] = 0;
      INDEXW[I] = 0;
    }
  }

// **
  // Loop over unreduced blocks
  IBEGIN = 1;
  WBEGIN = 1;
  for (JBLK = 1; JBLK <= NSPLIT.value; JBLK++) {
    IEND = ISPLIT[JBLK];
    IN = IEND - IBEGIN + 1;

    // 1 X 1 block
    if (IN == 1) {
      if ((IRANGE == ALLRNG) ||
          ((IRANGE == VALRNG) &&
              (D[IBEGIN] > VL.value) &&
              (D[IBEGIN] <= VU.value)) ||
          ((IRANGE == INDRNG) && (IBLOCK[WBEGIN] == JBLK))) {
        M.value++;
        W[M.value] = D[IBEGIN];
        WERR[M.value] = ZERO;
        // The gap for a single block doesn't matter for the later
        // algorithm and is assigned an arbitrary large value
        WGAP[M.value] = ZERO;
        IBLOCK[M.value] = JBLK;
        INDEXW[M.value] = 1;
        WBEGIN++;
      }
      // E[ IEND ] holds the shift for the initial RRR
      E[IEND] = ZERO;
      IBEGIN = IEND + 1;
      continue;
    }

    // Blocks of size larger than 1x1

    // E[ IEND ] will hold the shift for the initial RRR, for now set it =0
    E[IEND] = ZERO;

    // Find local outer bounds GL,GU for the block
    GL = D[IBEGIN];
    GU = D[IBEGIN];
    for (I = IBEGIN; I <= IEND; I++) {
      GL = min(GERS[2 * I - 1], GL);
      GU = max(GERS[2 * I], GU);
    }
    SPDIAM = GU - GL;

    if (!((IRANGE == ALLRNG) && (!FORCEB))) {
      // Count the number of eigenvalues in the current block.
      MB = 0;
      for (I = WBEGIN; I <= MM.value; I++) {
        if (IBLOCK[I] == JBLK) {
          MB++;
        } else {
          break;
        }
      }
      // }

      if (MB == 0) {
        // No eigenvalue in the current block lies in the desired range
        // E[ IEND ] holds the shift for the initial RRR
        E[IEND] = ZERO;
        IBEGIN = IEND + 1;
        continue;
      } else {
        // Decide whether dqds or bisection is more efficient
        USEDQD = ((MB > FAC * IN) && (!FORCEB));
        WEND = WBEGIN + MB - 1;
        // Calculate gaps for the current block
        // In later stages, when representations for individual
        // eigenvalues are different, we use SIGMA = E[ IEND ].
        SIGMA = ZERO;
        for (I = WBEGIN; I <= WEND - 1; I++) {
          WGAP[I] = max(ZERO, W[I + 1] - WERR[I + 1] - (W[I] + WERR[I]));
        }
        WGAP[WEND] = max(ZERO, VU.value - SIGMA - (W[WEND] + WERR[WEND]));
        // Find local index of the first and last desired evalue.
        INDL = INDEXW[WBEGIN];
        INDU = INDEXW[WEND];
      }
    }
    if (((IRANGE == ALLRNG) && (!FORCEB)) || USEDQD) {
      // Case of DQDS
      // Find approximations to the extremal eigenvalues of the block
      dlarrk(IN, 1, GL, GU, D(IBEGIN), E2(IBEGIN), PIVMIN.value, RTL, TMP, TMP1,
          IINFO);
      if (IINFO.value != 0) {
        INFO.value = -1;
        return;
      }
      ISLEFT = max(
          GL,
          TMP.value -
              TMP1.value -
              HNDRD * EPS * (TMP.value - TMP1.value).abs());
      dlarrk(IN, IN, GL, GU, D(IBEGIN), E2(IBEGIN), PIVMIN.value, RTL, TMP,
          TMP1, IINFO);
      if (IINFO.value != 0) {
        INFO.value = -1;
        return;
      }
      ISRGHT = min(
          GU,
          TMP.value +
              TMP1.value +
              HNDRD * EPS * (TMP.value + TMP1.value).abs());
      // Improve the estimate of the spectral diameter
      SPDIAM = ISRGHT - ISLEFT;
    } else {
      // Case of bisection
      // Find approximations to the wanted extremal eigenvalues
      ISLEFT = max(
          GL,
          W[WBEGIN] -
              WERR[WBEGIN] -
              HNDRD * EPS * (W[WBEGIN] - WERR[WBEGIN]).abs());
      ISRGHT = min(GU,
          W[WEND] + WERR[WEND] + HNDRD * EPS * (W[WEND] + WERR[WEND]).abs());
    }

    // Decide whether the base representation for the current block
    // L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
    // should be on the left or the right end of the current block.
    // The strategy is to shift to the end which is "more populated"
    // Furthermore, decide whether to use DQDS for the computation of
    // the eigenvalue approximations at the end of DLARRE or bisection.
    // dqds is chosen if all eigenvalues are desired or the number of
    // eigenvalues to be computed is large compared to the blocksize.
    if ((IRANGE == ALLRNG) && (!FORCEB)) {
      // If all the eigenvalues have to be computed, we use dqd
      USEDQD = true;
      // INDL is the local index of the first eigenvalue to compute
      INDL = 1;
      INDU = IN;
      // MB =  number of eigenvalues to compute
      MB = IN;
      WEND = WBEGIN + MB - 1;
      // Define 1/4 and 3/4 points of the spectrum
      S1 = ISLEFT + FOURTH * SPDIAM;
      S2 = ISRGHT - FOURTH * SPDIAM;
    } else {
      // DLARRD has computed IBLOCK and INDEXW for each eigenvalue
      // approximation.
      // choose sigma
      if (USEDQD) {
        S1 = ISLEFT + FOURTH * SPDIAM;
        S2 = ISRGHT - FOURTH * SPDIAM;
      } else {
        TMP.value = min(ISRGHT, VU.value) - max(ISLEFT, VL.value);
        S1 = max(ISLEFT, VL.value) + FOURTH * TMP.value;
        S2 = min(ISRGHT, VU.value) - FOURTH * TMP.value;
      }
    }

    // Compute the negcount at the 1/4 and 3/4 points
    if (MB > 1) {
      dlarrc('T', IN, S1, S2, D(IBEGIN), E(IBEGIN), PIVMIN.value, CNT, CNT1,
          CNT2, IINFO);
    }

    if (MB == 1) {
      SIGMA = GL;
      SGNDEF = ONE;
    } else if (CNT1.value - INDL >= INDU - CNT2.value) {
      if ((IRANGE == ALLRNG) && (!FORCEB)) {
        SIGMA = max(ISLEFT, GL);
      } else if (USEDQD) {
        // use Gerschgorin bound as shift to get pos def matrix
        // for dqds
        SIGMA = ISLEFT;
      } else {
        // use approximation of the first desired eigenvalue of the
        // block as shift
        SIGMA = max(ISLEFT, VL.value);
      }
      SGNDEF = ONE;
    } else {
      if ((IRANGE == ALLRNG) && (!FORCEB)) {
        SIGMA = min(ISRGHT, GU);
      } else if (USEDQD) {
        // use Gerschgorin bound as shift to get neg def matrix
        // for dqds
        SIGMA = ISRGHT;
      } else {
        // use approximation of the first desired eigenvalue of the
        // block as shift
        SIGMA = min(ISRGHT, VU.value);
      }
      SGNDEF = -ONE;
    }

    // An initial SIGMA has been chosen that will be used for computing
    // T - SIGMA I = L D L^T
    // Define the increment TAU of the shift in case the initial shift
    // needs to be refined to obtain a factorization with not too much
    // element growth.
    if (USEDQD) {
      // The initial SIGMA was to the outer end of the spectrum
      // the matrix is definite and we need not retreat.
      TAU = SPDIAM * EPS * N + TWO * PIVMIN.value;
      TAU = max(TAU, TWO * EPS * (SIGMA).abs());
    } else {
      if (MB > 1) {
        CLWDTH = W[WEND] + WERR[WEND] - W[WBEGIN] - WERR[WBEGIN];
        AVGAP = (CLWDTH / (WEND - WBEGIN).toDouble()).abs();
        if (SGNDEF == ONE) {
          TAU = HALF * max(WGAP[WBEGIN], AVGAP);
          TAU = max(TAU, WERR[WBEGIN]);
        } else {
          TAU = HALF * max(WGAP[WEND - 1], AVGAP);
          TAU = max(TAU, WERR[WEND]);
        }
      } else {
        TAU = WERR[WBEGIN];
      }
    }
    var baseFound = false;
    for (IDUM = 1; IDUM <= MAXTRY; IDUM++) {
      // Compute L D L^T factorization of tridiagonal matrix T - sigma I.
      // Store D in WORK[1:IN], L in WORK[IN+1:2*IN], and reciprocals of
      // pivots in WORK[2*IN+1:3*IN]
      DPIVOT = D[IBEGIN] - SIGMA;
      WORK[1] = DPIVOT;
      DMAX = (WORK[1]).abs();
      J = IBEGIN;
      for (I = 1; I <= IN - 1; I++) {
        WORK[2 * IN + I] = ONE / WORK[I];
        TMP.value = E[J] * WORK[2 * IN + I];
        WORK[IN + I] = TMP.value;
        DPIVOT = (D[J + 1] - SIGMA) - TMP.value * E[J];
        WORK[I + 1] = DPIVOT;
        DMAX = max(DMAX, (DPIVOT).abs());
        J++;
      }
      // check for element growth
      if (DMAX > MAXGROWTH * SPDIAM) {
        NOREP = true;
      } else {
        NOREP = false;
      }
      if (USEDQD && !NOREP) {
        // Ensure the definiteness of the representation
        // All entries of D (of L D L^T) must have the same sign
        for (I = 1; I <= IN; I++) {
          TMP.value = SGNDEF * WORK[I];
          if (TMP.value < ZERO) NOREP = true;
        }
      }
      if (NOREP) {
        // Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
        // shift which makes the matrix definite. So we should end up
        // here really only in the case of IRANGE = VALRNG or INDRNG.
        if (IDUM == MAXTRY - 1) {
          if (SGNDEF == ONE) {
            // The fudged Gerschgorin shift should succeed
            SIGMA = GL - FUDGE * SPDIAM * EPS * N - FUDGE * TWO * PIVMIN.value;
          } else {
            SIGMA = GU + FUDGE * SPDIAM * EPS * N + FUDGE * TWO * PIVMIN.value;
          }
        } else {
          SIGMA -= SGNDEF * TAU;
          TAU = TWO * TAU;
        }
      } else {
        // an initial RRR is found
        baseFound = true;
        break;
      }
    }
    if (!baseFound) {
      // if the program reaches this point, no base representation could be
      // found in MAXTRY iterations.
      INFO.value = 2;
      return;
    }

    // At this point, we have found an initial base representation
    // T - SIGMA I = L D L^T with not too much element growth.
    // Store the shift.
    E[IEND] = SIGMA;
    // Store D and L.
    dcopy(IN, WORK, 1, D(IBEGIN), 1);
    dcopy(IN - 1, WORK(IN + 1), 1, E(IBEGIN), 1);

    if (MB > 1) {
      // Perturb each entry of the base representation by a small
      // (but random) relative amount to overcome difficulties with
      // glued matrices.

      for (I = 1; I <= 4; I++) {
        ISEED[I] = 1;
      }

      dlarnv(2, ISEED, 2 * IN - 1, WORK(1));
      for (I = 1; I <= IN - 1; I++) {
        D[IBEGIN + I - 1] *= (ONE + EPS * PERT * WORK[I]);
        E[IBEGIN + I - 1] =
            E[IBEGIN + I - 1] * (ONE + EPS * PERT * WORK[IN + I]);
      }
      D[IEND] *= (ONE + EPS * FOUR * WORK[IN]);
    }

    // Don't update the Gerschgorin intervals because keeping track
    // of the updates would be too much work in DLARRV.
    // We update W instead and use it to locate the proper Gerschgorin
    // intervals.

    // Compute the required eigenvalues of L D L' by bisection or dqds
    if (!USEDQD) {
      // If DLARRD has been used, shift the eigenvalue approximations
      // according to their representation. This is necessary for
      // a uniform DLARRV since dqds computes eigenvalues of the
      // shifted representation. In DLARRV, W will always hold the
      // UNshifted eigenvalue approximation.
      for (J = WBEGIN; J <= WEND; J++) {
        W[J] -= SIGMA;
        WERR[J] += (W[J]).abs() * EPS;
      }
      // call DLARRB to reduce eigenvalue error of the approximations
      // from DLARRD
      for (I = IBEGIN; I <= IEND - 1; I++) {
        WORK[I] = D[I] * pow(E[I], 2);
      }
      // use bisection to find EV from INDL to INDU
      dlarrb(
          IN,
          D(IBEGIN),
          WORK(IBEGIN),
          INDL,
          INDU,
          RTOL1,
          RTOL2,
          INDL - 1,
          W(WBEGIN),
          WGAP(WBEGIN),
          WERR(WBEGIN),
          WORK(2 * N + 1),
          IWORK,
          PIVMIN.value,
          SPDIAM,
          IN,
          IINFO);
      if (IINFO.value != 0) {
        INFO.value = -4;
        return;
      }
      // DLARRB computes all gaps correctly except for the last one
      // Record distance to VU.value/GU
      WGAP[WEND] = max(ZERO, (VU.value - SIGMA) - (W[WEND] + WERR[WEND]));
      for (I = INDL; I <= INDU; I++) {
        M.value++;
        IBLOCK[M.value] = JBLK;
        INDEXW[M.value] = I;
      }
    } else {
      // Call dqds to get all eigs (and then possibly delete unwanted
      // eigenvalues).
      // Note that dqds finds the eigenvalues of the L D L^T representation
      // of T to high relative accuracy. High relative accuracy
      // might be lost when the shift of the RRR is subtracted to obtain
      // the eigenvalues of T. However, T is not guaranteed to define its
      // eigenvalues to high relative accuracy anyway.
      // Set RTOL to the order of the tolerance used in DLASQ2
      // This is an ESTIMATED error, the worst case bound is 4*N*EPS
      // which is usually too large and requires unnecessary work to be
      // done by bisection when computing the eigenvectors
      RTOL = log(IN) * FOUR * EPS;
      J = IBEGIN;
      for (I = 1; I <= IN - 1; I++) {
        WORK[2 * I - 1] = (D[J]).abs();
        WORK[2 * I] = E[J] * E[J] * WORK[2 * I - 1];
        J++;
      }
      WORK[2 * IN - 1] = (D[IEND]).abs();
      WORK[2 * IN] = ZERO;
      dlasq2(IN, WORK, IINFO);
      if (IINFO.value != 0) {
        // If IINFO = -5 then an index is part of a tight cluster
        // and should be changed. The index is in IWORK[1] and the
        // gap is in WORK[N+1]
        INFO.value = -5;
        return;
      } else {
        // Test that all eigenvalues are positive as expected
        for (I = 1; I <= IN; I++) {
          if (WORK[I] < ZERO) {
            INFO.value = -6;
            return;
          }
        }
      }
      if (SGNDEF > ZERO) {
        for (I = INDL; I <= INDU; I++) {
          M.value++;
          W[M.value] = WORK[IN - I + 1];
          IBLOCK[M.value] = JBLK;
          INDEXW[M.value] = I;
        }
      } else {
        for (I = INDL; I <= INDU; I++) {
          M.value++;
          W[M.value] = -WORK[I];
          IBLOCK[M.value] = JBLK;
          INDEXW[M.value] = I;
        }
      }

      for (I = M.value - MB + 1; I <= M.value; I++) {
        // the value of RTOL below should be the tolerance in DLASQ2
        WERR[I] = RTOL * (W[I]).abs();
      }
      for (I = M.value - MB + 1; I <= M.value - 1; I++) {
        // compute the right gap between the intervals
        WGAP[I] = max(ZERO, W[I + 1] - WERR[I + 1] - (W[I] + WERR[I]));
      }
      WGAP[M.value] =
          max(ZERO, (VU.value - SIGMA) - (W[M.value] + WERR[M.value]));
    }
    // proceed with next block
    IBEGIN = IEND + 1;
    WBEGIN = WEND + 1;
  }
}
