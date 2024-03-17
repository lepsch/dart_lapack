import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaebz.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlarrd(
  final String RANGE,
  final String ORDER,
  final int N,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final Array<double> GERS_,
  final double RELTOL,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> E2_,
  final double PIVMIN,
  final int NSPLIT,
  final Array<int> ISPLIT_,
  final Box<int> M,
  final Array<double> W_,
  final Array<double> WERR_,
  final Box<double> WL,
  final Box<double> WU,
  final Array<int> IBLOCK_,
  final Array<int> INDEXW_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final GERS = GERS_.having();
  final D = D_.having();
  final E = E_.having();
  final E2 = E2_.having();
  final ISPLIT = ISPLIT_.having();
  final W = W_.having();
  final WERR = WERR_.having();
  final IBLOCK = IBLOCK_.having();
  final INDEXW = INDEXW_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = ONE / TWO, FUDGE = TWO;
  const ALLRNG = 1, VALRNG = 2, INDRNG = 3;
  bool NCNVRG, TOOFEW;
  int I,
      IB,
      IBEGIN,
      IDISCL,
      IDISCU,
      IE,
      IEND,
      IN,
      IOFF,
      IRANGE,
      ITMAX,
      ITMP1,
      ITMP2,
      IW,
      IWOFF,
      J,
      JBLK,
      JDISC,
      JE,
      JEE,
      NB,
      NWL,
      NWU;
  double ATOLI,
      EPS,
      GL,
      GU,
      RTOLI,
      TMP1,
      TMP2,
      TNORM,
      UFLOW,
      WKILL,
      WLU = 0,
      WUL = 0;
  final IDUMMA = Array<int>(1);
  final IINFO = Box(0), IOUT = Box(0), IM = Box(0);

  INFO.value = 0;
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
  } else {
    IRANGE = 0;
  }

  // Check for Errors

  if (IRANGE <= 0) {
    INFO.value = -1;
  } else if (!(lsame(ORDER, 'B') || lsame(ORDER, 'E'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (IRANGE == VALRNG) {
    if (VL >= VU) INFO.value = -5;
  } else if (IRANGE == INDRNG && (IL < 1 || IL > max(1, N))) {
    INFO.value = -6;
  } else if (IRANGE == INDRNG && (IU < min(N, IL) || IU > N)) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    return;
  }

  // Initialize error flags
  NCNVRG = false;
  TOOFEW = false;

  // Simplification:
  if (IRANGE == INDRNG && IL == 1 && IU == N) IRANGE = 1;

  // Get machine constants
  EPS = dlamch('P');
  UFLOW = dlamch('U');

  // Special Case when N=1
  // Treat case of 1x1 matrix for quick return;
  if (N == 1) {
    if ((IRANGE == ALLRNG) ||
        ((IRANGE == VALRNG) && (D[1] > VL) && (D[1] <= VU)) ||
        ((IRANGE == INDRNG) && (IL == 1) && (IU == 1))) {
      M.value = 1;
      W[1] = D[1];
      // The computation error of the eigenvalue is zero
      WERR[1] = ZERO;
      IBLOCK[1] = 1;
      INDEXW[1] = 1;
    }
    return;
  }

  // NB is the minimum vector length for vector bisection, or 0
  // if only scalar is to be done.
  NB = ilaenv(1, 'DSTEBZ', ' ', N, -1, -1, -1);
  if (NB <= 1) NB = 0;

  // Find global spectral radius
  GL = D[1];
  GU = D[1];
  for (I = 1; I <= N; I++) {
    GL = min(GL, GERS[2 * I - 1]);
    GU = max(GU, GERS[2 * I]);
  }
  // Compute global Gerschgorin bounds and spectral diameter
  TNORM = max(GL.abs(), GU.abs());
  GL -= FUDGE * TNORM * EPS * N - FUDGE * TWO * PIVMIN;
  GU += FUDGE * TNORM * EPS * N + FUDGE * TWO * PIVMIN;
  // [JAN/28/2009] remove the line below since SPDIAM variable not use
  // SPDIAM = GU - GL
  // Input arguments for DLAEBZ:
  // The relative tolerance.  An interval (a,b] lies within
  // "relative tolerance" if  b-a < RELTOL*max(|a|,|b|),
  RTOLI = RELTOL;
  // Set the absolute tolerance for interval convergence to zero to force
  // interval convergence based on relative size of the interval.
  // This is dangerous because intervals might not converge when RELTOL is
  // small. But at least a very small number should be selected so that for
  // strongly graded matrices, the code can get relatively accurate
  // eigenvalues.
  ATOLI = FUDGE * TWO * UFLOW + FUDGE * TWO * PIVMIN;

  if (IRANGE == INDRNG) {
    // RANGE='I': Compute an interval containing eigenvalues
    // IL through IU. The initial interval [GL,GU] from the global
    // Gerschgorin bounds GL and GU is refined by DLAEBZ.
    ITMAX = (log(TNORM + PIVMIN) - log(PIVMIN)) ~/ log(TWO) + 2;
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

    dlaebz(
        3,
        ITMAX,
        N,
        2,
        2,
        NB,
        ATOLI,
        RTOLI,
        PIVMIN,
        D,
        E,
        E2,
        IWORK(5),
        WORK(N + 1).asMatrix(2),
        WORK(N + 5),
        IOUT,
        IWORK.asMatrix(2),
        W,
        IBLOCK,
        IINFO);
    if (IINFO.value != 0) {
      INFO.value = IINFO.value;
      return;
    }
    // On exit, output intervals may not be ordered by ascending negcount
    if (IWORK[6] == IU) {
      WL.value = WORK[N + 1];
      WLU = WORK[N + 3];
      NWL = IWORK[1];
      WU.value = WORK[N + 4];
      WUL = WORK[N + 2];
      NWU = IWORK[4];
    } else {
      WL.value = WORK[N + 2];
      WLU = WORK[N + 4];
      NWL = IWORK[2];
      WU.value = WORK[N + 3];
      WUL = WORK[N + 1];
      NWU = IWORK[3];
    }
    // On exit, the interval [WL.value, WLU] contains a value with negcount NWL,
    // and [WUL, WU.value] contains a value with negcount NWU.
    if (NWL < 0 || NWL >= N || NWU < 1 || NWU > N) {
      INFO.value = 4;
      return;
    }
  } else if (IRANGE == VALRNG) {
    WL.value = VL;
    WU.value = VU;
  } else if (IRANGE == ALLRNG) {
    WL.value = GL;
    WU.value = GU;
  }

  // Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
  // NWL accumulates the number of eigenvalues <= WL.value,
  // NWU accumulates the number of eigenvalues <= WU.value
  M.value = 0;
  IEND = 0;
  INFO.value = 0;
  NWL = 0;
  NWU = 0;

  for (JBLK = 1; JBLK <= NSPLIT; JBLK++) {
    IOFF = IEND;
    IBEGIN = IOFF + 1;
    IEND = ISPLIT[JBLK];
    IN = IEND - IOFF;

    if (IN == 1) {
      // 1x1 block
      if (WL.value >= D[IBEGIN] - PIVMIN) NWL++;
      if (WU.value >= D[IBEGIN] - PIVMIN) NWU++;
      if (IRANGE == ALLRNG ||
          (WL.value < D[IBEGIN] - PIVMIN && WU.value >= D[IBEGIN] - PIVMIN)) {
        M.value++;
        W[M.value] = D[IBEGIN];
        WERR[M.value] = ZERO;
        // The gap for a single block doesn't matter for the later
        // algorithm and is assigned an arbitrary large value
        IBLOCK[M.value] = JBLK;
        INDEXW[M.value] = 1;
      }

      // Disabled 2x2 case because of a failure on the following matrix
      // RANGE = 'I', IL = IU = 4
      //   Original Tridiagonal, d = [
      //    -0.150102010615740e+00
      //    -0.849897989384260e+00
      //    -0.128208148052635e-15
      //     0.128257718286320e-15
      //   ];
      //   e = [
      //    -0.357171383266986e+00
      //    -0.180411241501588e-15
      //    -0.175152352710251e-15
      //   ];

      // ELSE if( IN == 2 ) THEN
      //    // 2x2 block
      //    DISC = sqrt( (HALF*(D[IBEGIN]-D[IEND]))**2 + E[IBEGIN]**2 )
      //    TMP1 = HALF*(D[IBEGIN]+D[IEND])
      //    L1 = TMP1 - DISC
      //    if( WL.value >= L1-PIVMIN )
      //        NWL++
      //    if( WU.value >= L1-PIVMIN )
      //        NWU++
      //    if( IRANGE == ALLRNG || ( WL.value < L1-PIVMIN && WU.value.GE.
      //        L1-PIVMIN ) ) THEN
      //        M.value++
      //        W[ M.value ] = L1
      //        // The uncertainty of eigenvalues of a 2x2 matrix is very small
      //        WERR[ M.value ] = EPS * ABS( W[ M.value ] ) * TWO
      //        IBLOCK[ M.value ] = JBLK
      //        INDEXW[ M.value ] = 1
      //     ENDIF
      //     L2 = TMP1 + DISC
      //     if( WL.value >= L2-PIVMIN )
      //         NWL++
      //     if( WU.value >= L2-PIVMIN )
      //         NWU++
      //     if( IRANGE == ALLRNG || ( WL.value < L2-PIVMIN && WU.value.GE.
      //         L2-PIVMIN ) ) THEN
      //        M.value++
      //        W[ M.value ] = L2
      //        // The uncertainty of eigenvalues of a 2x2 matrix is very small
      //        WERR[ M.value ] = EPS * ABS( W[ M.value ] ) * TWO
      //        IBLOCK[ M.value ] = JBLK
      //        INDEXW[ M.value ] = 2
      //     ENDIF
    } else {
      // General Case - block of size IN >= 2
      // Compute local Gerschgorin interval and use it as the initial
      // interval for DLAEBZ
      GU = D[IBEGIN];
      GL = D[IBEGIN];
      TMP1 = ZERO;

      for (J = IBEGIN; J <= IEND; J++) {
        GL = min(GL, GERS[2 * J - 1]);
        GU = max(GU, GERS[2 * J]);
      }
      // [JAN/28/2009]
      // change SPDIAM by TNORM in lines 2 and 3 thereafter
      // line 1: remove computation of SPDIAM (not useful anymore)
      // SPDIAM = GU - GL
      // GL -= FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN
      // GU += FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN
      GL -= FUDGE * TNORM * EPS * IN - FUDGE * PIVMIN;
      GU += FUDGE * TNORM * EPS * IN + FUDGE * PIVMIN;

      if (IRANGE > 1) {
        if (GU < WL.value) {
          // the local block contains none of the wanted eigenvalues
          NWL += IN;
          NWU += IN;
          continue;
        }
        // refine search interval if possible, only range (WL.value,WU.value] matters
        GL = max(GL, WL.value);
        GU = min(GU, WU.value);
        if (GL >= GU) continue;
      }

      // Find negcount of initial interval boundaries GL and GU
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
          D(IBEGIN),
          E(IBEGIN),
          E2(IBEGIN),
          IDUMMA,
          WORK(N + 1).asMatrix(IN),
          WORK(N + 2 * IN + 1),
          IM,
          IWORK.asMatrix(IN),
          W(M.value + 1),
          IBLOCK(M.value + 1),
          IINFO);
      if (IINFO.value != 0) {
        INFO.value = IINFO.value;
        return;
      }

      NWL += IWORK[1];
      NWU += IWORK[IN + 1];
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
          D(IBEGIN),
          E(IBEGIN),
          E2(IBEGIN),
          IDUMMA,
          WORK(N + 1).asMatrix(IN),
          WORK(N + 2 * IN + 1),
          IOUT,
          IWORK.asMatrix(IN),
          W(M.value + 1),
          IBLOCK(M.value + 1),
          IINFO);
      if (IINFO.value != 0) {
        INFO.value = IINFO.value;
        return;
      }

      // Copy eigenvalues into W and IBLOCK
      // Use -JBLK for block number for unconverged eigenvalues.
      // Loop over the number of output intervals from DLAEBZ
      for (J = 1; J <= IOUT.value; J++) {
        // eigenvalue approximation is middle point of interval
        TMP1 = HALF * (WORK[J + N] + WORK[J + IN + N]);
        // semi length of error interval
        TMP2 = HALF * (WORK[J + N] - WORK[J + IN + N]).abs();
        if (J > IOUT.value - IINFO.value) {
          // Flag non-convergence.
          NCNVRG = true;
          IB = -JBLK;
        } else {
          IB = JBLK;
        }
        for (JE = IWORK[J] + 1 + IWOFF; JE <= IWORK[J + IN] + IWOFF; JE++) {
          W[JE] = TMP1;
          WERR[JE] = TMP2;
          INDEXW[JE] = JE - IWOFF;
          IBLOCK[JE] = IB;
        }
      }

      M.value += IM.value;
    }
  }

  // If RANGE='I', then (WL.value,WU.value) contains eigenvalues NWL+1,...,NWU
  // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
  if (IRANGE == INDRNG) {
    IDISCL = IL - 1 - NWL;
    IDISCU = NWU - IU;

    if (IDISCL > 0) {
      IM.value = 0;
      for (JE = 1; JE <= M.value; JE++) {
        // Remove some of the smallest eigenvalues from the left so that
        // at the end IDISCL =0. Move all eigenvalues up to the left.
        if (W[JE] <= WLU && IDISCL > 0) {
          IDISCL--;
        } else {
          IM.value++;
          W[IM.value] = W[JE];
          WERR[IM.value] = WERR[JE];
          INDEXW[IM.value] = INDEXW[JE];
          IBLOCK[IM.value] = IBLOCK[JE];
        }
      }
      M.value = IM.value;
    }
    if (IDISCU > 0) {
      // Remove some of the largest eigenvalues from the right so that
      // at the end IDISCU =0. Move all eigenvalues up to the left.
      IM.value = M.value + 1;
      for (JE = M.value; JE >= 1; JE--) {
        if (W[JE] >= WUL && IDISCU > 0) {
          IDISCU--;
        } else {
          IM.value--;
          W[IM.value] = W[JE];
          WERR[IM.value] = WERR[JE];
          INDEXW[IM.value] = INDEXW[JE];
          IBLOCK[IM.value] = IBLOCK[JE];
        }
      }
      JEE = 0;
      for (JE = IM.value; JE <= M.value; JE++) {
        JEE++;
        W[JEE] = W[JE];
        WERR[JEE] = WERR[JE];
        INDEXW[JEE] = INDEXW[JE];
        IBLOCK[JEE] = IBLOCK[JE];
      }
      M.value -= IM.value + 1;
    }

    if (IDISCL > 0 || IDISCU > 0) {
      // Code to deal with effects of bad arithmetic. (If N(w) is
      // monotone non-decreasing, this should never happen.)
      // Some low eigenvalues to be discarded are not in (WL.value,WLU],
      // or high eigenvalues to be discarded are not in (WUL,WU.value]
      // so just kill off the smallest IDISCL/largest IDISCU
      // eigenvalues, by marking the corresponding IBLOCK = 0
      if (IDISCL > 0) {
        WKILL = WU.value;
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
        WKILL = WL.value;
        for (JDISC = 1; JDISC <= IDISCU; JDISC++) {
          IW = 0;
          for (JE = 1; JE <= M.value; JE++) {
            if (IBLOCK[JE] != 0 && (W[JE] >= WKILL || IW == 0)) {
              IW = JE;
              WKILL = W[JE];
            }
          }
          IBLOCK[IW] = 0;
        }
      }
      // Now erase all eigenvalues with IBLOCK set to zero
      IM.value = 0;
      for (JE = 1; JE <= M.value; JE++) {
        if (IBLOCK[JE] != 0) {
          IM.value++;
          W[IM.value] = W[JE];
          WERR[IM.value] = WERR[JE];
          INDEXW[IM.value] = INDEXW[JE];
          IBLOCK[IM.value] = IBLOCK[JE];
        }
      }
      M.value = IM.value;
    }
    if (IDISCL < 0 || IDISCU < 0) {
      TOOFEW = true;
    }
  }

  if ((IRANGE == ALLRNG && M.value != N) ||
      (IRANGE == INDRNG && M.value != IU - IL + 1)) {
    TOOFEW = true;
  }

  // If ORDER='B', do nothing the eigenvalues are already sorted by
  // block.
  // If ORDER='E', sort the eigenvalues from smallest to largest

  if (lsame(ORDER, 'E') && NSPLIT > 1) {
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
        TMP2 = WERR[IE];
        ITMP1 = IBLOCK[IE];
        ITMP2 = INDEXW[IE];
        W[IE] = W[JE];
        WERR[IE] = WERR[JE];
        IBLOCK[IE] = IBLOCK[JE];
        INDEXW[IE] = INDEXW[JE];
        W[JE] = TMP1;
        WERR[JE] = TMP2;
        IBLOCK[JE] = ITMP1;
        INDEXW[JE] = ITMP2;
      }
    }
  }

  INFO.value = 0;
  if (NCNVRG) INFO.value++;
  if (TOOFEW) INFO.value += 2;
}
