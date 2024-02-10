import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaebz.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dstebz(
  final String RANGE,
  final String ORDER,
  final int N,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Array<double> D,
  final Array<double> E,
  final Box<int> M,
  final Box<int> NSPLIT,
  final Array<double> W,
  final Array<int> IBLOCK,
  final Array<int> ISPLIT,
  final Array<double> WORK,
  final Array<int> IWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 1.0 / TWO;
  const FUDGE = 2.1, RELFAC = 2.0;
  bool NCNVRG, TOOFEW;
  int IB,
      IBEGIN,
      IDISCL,
      IDISCU,
      IE,
      IEND,
      IM = 0,
      IN,
      IOFF,
      IORDER,
      IOUT = 0,
      IRANGE,
      ITMAX,
      ITMP1,
      IW,
      IWOFF,
      J,
      JB,
      JDISC,
      JE,
      NB,
      NWL,
      NWU;
  double ATOLI,
      BNORM,
      GL,
      GU,
      PIVMIN,
      RTOLI,
      SAFEMN,
      TMP1,
      TMP2,
      TNORM,
      ULP,
      WKILL,
      WL,
      WLU = 0,
      WU = 0,
      WUL = 0;
  final IDUMMA = Array<int>(1);
  final IINFO = Box(0);

  INFO.value = 0;

  // Decode RANGE

  if (lsame(RANGE, 'A')) {
    IRANGE = 1;
  } else if (lsame(RANGE, 'V')) {
    IRANGE = 2;
  } else if (lsame(RANGE, 'I')) {
    IRANGE = 3;
  } else {
    IRANGE = 0;
  }

  // Decode ORDER

  if (lsame(ORDER, 'B')) {
    IORDER = 2;
  } else if (lsame(ORDER, 'E')) {
    IORDER = 1;
  } else {
    IORDER = 0;
  }

  // Check for Errors

  if (IRANGE <= 0) {
    INFO.value = -1;
  } else if (IORDER <= 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (IRANGE == 2) {
    if (VL >= VU) INFO.value = -5;
  } else if (IRANGE == 3 && (IL < 1 || IL > max(1, N))) {
    INFO.value = -6;
  } else if (IRANGE == 3 && (IU < min(N, IL) || IU > N)) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    xerbla('DSTEBZ', -INFO.value);
    return;
  }

  // Initialize error flags

  INFO.value = 0;
  NCNVRG = false;
  TOOFEW = false;

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  // Simplifications:

  if (IRANGE == 3 && IL == 1 && IU == N) IRANGE = 1;

  // Get machine constants
  // NB is the minimum vector length for vector bisection, or 0
  // if only scalar is to be done.

  SAFEMN = dlamch('S');
  ULP = dlamch('P');
  RTOLI = ULP * RELFAC;
  NB = ilaenv(1, 'DSTEBZ', ' ', N, -1, -1, -1);
  if (NB <= 1) NB = 0;

  // Special Case when N=1

  if (N == 1) {
    NSPLIT.value = 1;
    ISPLIT[1] = 1;
    if (IRANGE == 2 && (VL >= D[1] || VU < D[1])) {
      M.value = 0;
    } else {
      W[1] = D[1];
      IBLOCK[1] = 1;
      M.value = 1;
    }
    return;
  }

  // Compute Splitting Points

  NSPLIT.value = 1;
  WORK[N] = ZERO;
  PIVMIN = ONE;

  for (J = 2; J <= N; J++) {
    TMP1 = pow(E[J - 1], 2).toDouble();
    if ((D[J] * D[J - 1]).abs() * pow(ULP, 2) + SAFEMN > TMP1) {
      ISPLIT[NSPLIT.value] = J - 1;
      NSPLIT.value = NSPLIT.value + 1;
      WORK[J - 1] = ZERO;
    } else {
      WORK[J - 1] = TMP1;
      PIVMIN = max(PIVMIN, TMP1);
    }
  }
  ISPLIT[NSPLIT.value] = N;
  PIVMIN = PIVMIN * SAFEMN;

  // Compute Interval and ATOLI

  if (IRANGE == 3) {
    // RANGE='I': Compute the interval containing eigenvalues
    // IL through IU.

    // Compute Gershgorin interval for entire (split) matrix
    // and use it as the initial interval

    GU = D[1];
    GL = D[1];
    TMP1 = ZERO;

    for (J = 1; J <= N - 1; J++) {
      TMP2 = sqrt(WORK[J]);
      GU = max(GU, D[J] + TMP1 + TMP2);
      GL = min(GL, D[J] - TMP1 - TMP2);
      TMP1 = TMP2;
    }

    GU = max(GU, D[N] + TMP1);
    GL = min(GL, D[N] - TMP1);
    TNORM = max((GL).abs(), (GU).abs());
    GL = GL - FUDGE * TNORM * ULP * N - FUDGE * TWO * PIVMIN;
    GU = GU + FUDGE * TNORM * ULP * N + FUDGE * PIVMIN;

    // Compute Iteration parameters

    ITMAX = (log(TNORM + PIVMIN) - log(PIVMIN)) ~/ log(TWO) + 2;
    if (ABSTOL <= ZERO) {
      ATOLI = ULP * TNORM;
    } else {
      ATOLI = ABSTOL;
    }

    WORK[N + 1] = GL;
    WORK[N + 2] = GL;
    WORK[N + 3] = GU;
    WORK[N + 4] = GU;
    WORK[N + 5] = GL;
    WORK[N + 6] = GU;
    IWORK[1] = -1;
    IWORK[2] = -1;
    IWORK[3] = N + 1;
    IWORK[4] = N + 1;
    IWORK[5] = IL - 1;
    IWORK[6] = IU;

    dlaebz(3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, WORK, IWORK[5],
        WORK[N + 1], WORK[N + 5], IOUT, IWORK, W, IBLOCK, IINFO.value);

    if (IWORK[6] == IU) {
      WL = WORK[N + 1];
      WLU = WORK[N + 3];
      NWL = IWORK[1];
      WU = WORK[N + 4];
      WUL = WORK[N + 2];
      NWU = IWORK[4];
    } else {
      WL = WORK[N + 2];
      WLU = WORK[N + 4];
      NWL = IWORK[2];
      WU = WORK[N + 3];
      WUL = WORK[N + 1];
      NWU = IWORK[3];
    }

    if (NWL < 0 || NWL >= N || NWU < 1 || NWU > N) {
      INFO.value = 4;
      return;
    }
  } else {
    // RANGE='A' or 'V' -- Set ATOLI

    TNORM = max((D[1]).abs() + (E[1]).abs(), (D[N]).abs() + (E[N - 1]).abs());

    for (J = 2; J <= N - 1; J++) {
      TNORM = max(TNORM, (D[J]).abs() + (E[J - 1]).abs() + (E[J]).abs());
    }

    if (ABSTOL <= ZERO) {
      ATOLI = ULP * TNORM;
    } else {
      ATOLI = ABSTOL;
    }

    if (IRANGE == 2) {
      WL = VL;
      WU = VU;
    } else {
      WL = ZERO;
      WU = ZERO;
    }
  }

  // Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
  // NWL accumulates the number of eigenvalues <= WL,
  // NWU accumulates the number of eigenvalues <= WU

  M.value = 0;
  IEND = 0;
  INFO.value = 0;
  NWL = 0;
  NWU = 0;

  for (JB = 1; JB <= NSPLIT.value; JB++) {
    IOFF = IEND;
    IBEGIN = IOFF + 1;
    IEND = ISPLIT[JB];
    IN = IEND - IOFF;

    if (IN == 1) {
      // Special Case -- IN=1

      if (IRANGE == 1 || WL >= D[IBEGIN] - PIVMIN) NWL = NWL + 1;
      if (IRANGE == 1 || WU >= D[IBEGIN] - PIVMIN) NWU = NWU + 1;
      if (IRANGE == 1 ||
          (WL < D[IBEGIN] - PIVMIN && WU >= D[IBEGIN] - PIVMIN)) {
        M.value = M.value + 1;
        W[M.value] = D[IBEGIN];
        IBLOCK[M.value] = JB;
      }
    } else {
      // General Case -- IN > 1

      // Compute Gershgorin Interval
      // and use it as the initial interval

      GU = D[IBEGIN];
      GL = D[IBEGIN];
      TMP1 = ZERO;

      for (J = IBEGIN; J <= IEND - 1; J++) {
        TMP2 = (E[J]).abs();
        GU = max(GU, D[J] + TMP1 + TMP2);
        GL = min(GL, D[J] - TMP1 - TMP2);
        TMP1 = TMP2;
      }

      GU = max(GU, D[IEND] + TMP1);
      GL = min(GL, D[IEND] - TMP1);
      BNORM = max((GL).abs(), (GU).abs());
      GL = GL - FUDGE * BNORM * ULP * IN - FUDGE * PIVMIN;
      GU = GU + FUDGE * BNORM * ULP * IN + FUDGE * PIVMIN;

      // Compute ATOLI for the current submatrix

      if (ABSTOL <= ZERO) {
        ATOLI = ULP * max((GL).abs(), (GU).abs());
      } else {
        ATOLI = ABSTOL;
      }

      if (IRANGE > 1) {
        if (GU < WL) {
          NWL = NWL + IN;
          NWU = NWU + IN;
          continue;
        }
        GL = max(GL, WL);
        GU = min(GU, WU);
        if (GL >= GU) continue;
      }

      // Set Up Initial Interval

      WORK[N + 1] = GL;
      WORK[N + IN + 1] = GU;
      dlaebz(
          1,
          0,
          IN,
          IN,
          1,
          NB,
          ATOLI,
          RTOLI,
          PIVMIN,
          D[IBEGIN],
          E[IBEGIN],
          WORK[IBEGIN],
          IDUMMA,
          WORK[N + 1],
          WORK[N + 2 * IN + 1],
          IM,
          IWORK,
          W[M.value + 1],
          IBLOCK[M.value + 1],
          IINFO.value);

      NWL = NWL + IWORK[1];
      NWU = NWU + IWORK[IN + 1];
      IWOFF = M.value - IWORK[1];

      // Compute Eigenvalues

      ITMAX = (log(GU - GL + PIVMIN) - log(PIVMIN)) ~/ log(TWO) + 2;
      dlaebz(
          2,
          ITMAX,
          IN,
          IN,
          1,
          NB,
          ATOLI,
          RTOLI,
          PIVMIN,
          D[IBEGIN],
          E[IBEGIN],
          WORK[IBEGIN],
          IDUMMA,
          WORK[N + 1],
          WORK[N + 2 * IN + 1],
          IOUT,
          IWORK,
          W[M.value + 1],
          IBLOCK[M.value + 1],
          IINFO.value);

      // Copy Eigenvalues Into W and IBLOCK
      // Use -JB for block number for unconverged eigenvalues.

      for (J = 1; J <= IOUT; J++) {
        TMP1 = HALF * (WORK[J + N] + WORK[J + IN + N]);

        // Flag non-convergence.

        if (J > IOUT - IINFO.value) {
          NCNVRG = true;
          IB = -JB;
        } else {
          IB = JB;
        }
        for (JE = IWORK[J] + 1 + IWOFF; JE <= IWORK[J + IN] + IWOFF; JE++) {
          W[JE] = TMP1;
          IBLOCK[JE] = IB;
        }
      }

      M.value = M.value + IM;
    }
  }

  // If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
  // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.

  if (IRANGE == 3) {
    IM = 0;
    IDISCL = IL - 1 - NWL;
    IDISCU = NWU - IU;

    if (IDISCL > 0 || IDISCU > 0) {
      for (JE = 1; JE <= M.value; JE++) {
        if (W[JE] <= WLU && IDISCL > 0) {
          IDISCL = IDISCL - 1;
        } else if (W[JE] >= WUL && IDISCU > 0) {
          IDISCU = IDISCU - 1;
        } else {
          IM = IM + 1;
          W[IM] = W[JE];
          IBLOCK[IM] = IBLOCK[JE];
        }
      }
      M.value = IM;
    }
    if (IDISCL > 0 || IDISCU > 0) {
      // Code to deal with effects of bad arithmetic:
      // Some low eigenvalues to be discarded are not in (WL,WLU],
      // or high eigenvalues to be discarded are not in (WUL,WU]
      // so just kill off the smallest IDISCL/largest IDISCU
      // eigenvalues, by simply finding the smallest/largest
      // eigenvalue(s).

      // (If N(w) is monotone non-decreasing, this should never
      // happen.)

      if (IDISCL > 0) {
        WKILL = WU;
        for (JDISC = 1; JDISC <= IDISCL; JDISC++) {
          IW = 0;
          for (JE = 1; JE <= M.value; JE++) {
            if (IBLOCK[JE] != 0 && (W[JE] < WKILL || IW == 0)) {
              IW = JE;
              WKILL = W[JE];
            }
          }
          IBLOCK[IW] = 0;
        }
      }
      if (IDISCU > 0) {
        WKILL = WL;
        for (JDISC = 1; JDISC <= IDISCU; JDISC++) {
          IW = 0;
          for (JE = 1; JE <= M.value; JE++) {
            if (IBLOCK[JE] != 0 && (W[JE] > WKILL || IW == 0)) {
              IW = JE;
              WKILL = W[JE];
            }
          }
          IBLOCK[IW] = 0;
        }
      }
      IM = 0;
      for (JE = 1; JE <= M.value; JE++) {
        if (IBLOCK[JE] != 0) {
          IM = IM + 1;
          W[IM] = W[JE];
          IBLOCK[IM] = IBLOCK[JE];
        }
      }
      M.value = IM;
    }
    if (IDISCL < 0 || IDISCU < 0) {
      TOOFEW = true;
    }
  }

  // If ORDER='B', do nothing -- the eigenvalues are already sorted
  // by block.
  // If ORDER='E', sort the eigenvalues from smallest to largest

  if (IORDER == 1 && NSPLIT.value > 1) {
    for (JE = 1; JE <= M.value - 1; JE++) {
      IE = 0;
      TMP1 = W[JE];
      for (J = JE + 1; J <= M.value; J++) {
        if (W[J] < TMP1) {
          IE = J;
          TMP1 = W[J];
        }
      }

      if (IE != 0) {
        ITMP1 = IBLOCK[IE];
        W[IE] = W[JE];
        IBLOCK[IE] = IBLOCK[JE];
        W[JE] = TMP1;
        IBLOCK[JE] = ITMP1;
      }
    }
  }

  INFO.value = 0;
  if (NCNVRG) INFO.value = INFO.value + 1;
  if (TOOFEW) INFO.value = INFO.value + 2;
}
