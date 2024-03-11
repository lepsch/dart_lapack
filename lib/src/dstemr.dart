import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlae2.dart';
import 'package:lapack/src/dlaev2.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlarrc.dart';
import 'package:lapack/src/dlarre.dart';
import 'package:lapack/src/dlarrj.dart';
import 'package:lapack/src/dlarrr.dart';
import 'package:lapack/src/dlarrv.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dstemr(
  final String JOBZ,
  final String RANGE,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final int NZC,
  final Array<int> ISUPPZ_,
  final Box<bool> TRYRAC,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final ISUPPZ = ISUPPZ_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, FOUR = 4.0, MINRGP = 1.0e-3;
  bool ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ, ZQUERY, LAESWAP;
  int I,
      IBEGIN,
      IEND,
      IFIRST,
      IIL,
      IINDBL,
      IINDW,
      IINDWK,
      IINSPL,
      IIU,
      ILAST,
      IN,
      INDD,
      INDE2,
      INDERR,
      INDGP,
      INDGRS,
      INDWRK,
      J,
      JBLK,
      JJ,
      LIWMIN,
      LWMIN,
      OFFSET,
      WBEGIN,
      WEND = 0;
  double BIGNUM, EPS, RMAX, RMIN, SAFMIN, SCALE, SMLNUM, THRESH, TMP, TNRM;
  final WL = Box(0.0),
      WU = Box(0.0),
      PIVMIN = Box(0.0),
      RTOL1 = Box(0.0),
      RTOL2 = Box(0.0),
      CS = Box(0.0),
      SN = Box(0.0),
      R1 = Box(0.0),
      R2 = Box(0.0);

  final IINFO = Box(0),
      NSPLIT = Box(0),
      ITMP = Box(0),
      ITMP2 = Box(0),
      NZCMIN = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

  LQUERY = ((LWORK == -1) || (LIWORK == -1));
  ZQUERY = (NZC == -1);
  LAESWAP = false;

  // DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
  // In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
  // Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N.
  if (WANTZ) {
    LWMIN = 18 * N;
    LIWMIN = 10 * N;
  } else {
    // need less workspace if only the eigenvalues are wanted
    LWMIN = 12 * N;
    LIWMIN = 8 * N;
  }

  WL.value = ZERO;
  WU.value = ZERO;
  IIL = 0;
  IIU = 0;
  NSPLIT.value = 0;

  if (VALEIG) {
    // We do not reference VL, VU in the cases RANGE = 'I','A'
    // The interval (WL.value, WU.value] contains all the wanted eigenvalues.
    // It is either given by the user or computed in DLARRE.
    WL.value = VL;
    WU.value = VU;
  } else if (INDEIG) {
    // We do not reference IL, IU in the cases RANGE = 'V','A'
    IIL = IL;
    IIU = IU;
  }

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (VALEIG && N > 0 && WU.value <= WL.value) {
    INFO.value = -7;
  } else if (INDEIG && (IIL < 1 || IIL > N)) {
    INFO.value = -8;
  } else if (INDEIG && (IIU < IIL || IIU > N)) {
    INFO.value = -9;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -13;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -17;
  } else if (LIWORK < LIWMIN && !LQUERY) {
    INFO.value = -19;
  }

  // Get machine constants.

  SAFMIN = dlamch('Safe minimum');
  EPS = dlamch('Precision');
  SMLNUM = SAFMIN / EPS;
  BIGNUM = ONE / SMLNUM;
  RMIN = sqrt(SMLNUM);
  RMAX = min(sqrt(BIGNUM), ONE / sqrt(sqrt(SAFMIN)));

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (WANTZ && ALLEIG) {
      NZCMIN.value = N;
    } else if (WANTZ && VALEIG) {
      dlarrc('T', N, VL, VU, D, E, SAFMIN, NZCMIN, ITMP, ITMP2, INFO);
    } else if (WANTZ && INDEIG) {
      NZCMIN.value = IIU - IIL + 1;
    } else {
      // WANTZ == FALSE.
      NZCMIN.value = 0;
    }
    if (ZQUERY && INFO.value == 0) {
      Z[1][1] = NZCMIN.value.toDouble();
    } else if (NZC < NZCMIN.value && !ZQUERY) {
      INFO.value = -14;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSTEMR', -INFO.value);

    return;
  } else if (LQUERY || ZQUERY) {
    return;
  }

  // Handle N = 0, 1, and 2 cases immediately

  M.value = 0;
  if (N == 0) return;

  if (N == 1) {
    if (ALLEIG || INDEIG) {
      M.value = 1;
      W[1] = D[1];
    } else {
      if (WL.value < D[1] && WU.value >= D[1]) {
        M.value = 1;
        W[1] = D[1];
      }
    }
    if (WANTZ && (!ZQUERY)) {
      Z[1][1] = ONE;
      ISUPPZ[1] = 1;
      ISUPPZ[2] = 1;
    }
    return;
  }

  if (N == 2) {
    if (!WANTZ) {
      dlae2(D[1], E[1], D[2], R1, R2);
    } else if (WANTZ && (!ZQUERY)) {
      dlaev2(D[1], E[1], D[2], R1, R2, CS, SN);
    }
    // D/S/LAE2 and D/S/LAEV2 outputs satisfy |R1.value| >= |R2.value|. However,
    // the following code requires R1.value >= R2.value. Hence, we correct
    // the order of R1.value, R2.value, CS, SN if R1.value < R2.value before further processing.
    if (R1.value < R2.value) {
      E[2] = R1.value;
      R1.value = R2.value;
      R2.value = E[2];
      LAESWAP = true;
    }
    if (ALLEIG ||
        (VALEIG && (R2.value > WL.value) && (R2.value <= WU.value)) ||
        (INDEIG && (IIL == 1))) {
      M.value = M.value + 1;
      W[M.value] = R2.value;
      if (WANTZ && (!ZQUERY)) {
        if (LAESWAP) {
          Z[1][M.value] = CS.value;
          Z[2][M.value] = SN.value;
        } else {
          Z[1][M.value] = -SN.value;
          Z[2][M.value] = CS.value;
        }
        // Note: At most one of SN and CS can be zero.
        if (SN.value != ZERO) {
          if (CS.value != ZERO) {
            ISUPPZ[2 * M.value - 1] = 1;
            ISUPPZ[2 * M.value] = 2;
          } else {
            ISUPPZ[2 * M.value - 1] = 1;
            ISUPPZ[2 * M.value] = 1;
          }
        } else {
          ISUPPZ[2 * M.value - 1] = 2;
          ISUPPZ[2 * M.value] = 2;
        }
      }
    }
    if (ALLEIG ||
        (VALEIG && (R1.value > WL.value) && (R1.value <= WU.value)) ||
        (INDEIG && (IIU == 2))) {
      M.value = M.value + 1;
      W[M.value] = R1.value;
      if (WANTZ && (!ZQUERY)) {
        if (LAESWAP) {
          Z[1][M.value] = -SN.value;
          Z[2][M.value] = CS.value;
        } else {
          Z[1][M.value] = CS.value;
          Z[2][M.value] = SN.value;
        }
        // Note: At most one of SN and CS can be zero.
        if (SN.value != ZERO) {
          if (CS.value != ZERO) {
            ISUPPZ[2 * M.value - 1] = 1;
            ISUPPZ[2 * M.value] = 2;
          } else {
            ISUPPZ[2 * M.value - 1] = 1;
            ISUPPZ[2 * M.value] = 1;
          }
        } else {
          ISUPPZ[2 * M.value - 1] = 2;
          ISUPPZ[2 * M.value] = 2;
        }
      }
    }
  } else {
    // Continue with general N

    INDGRS = 1;
    INDERR = 2 * N + 1;
    INDGP = 3 * N + 1;
    INDD = 4 * N + 1;
    INDE2 = 5 * N + 1;
    INDWRK = 6 * N + 1;

    IINSPL = 1;
    IINDBL = N + 1;
    IINDW = 2 * N + 1;
    IINDWK = 3 * N + 1;

    // Scale matrix to allowable range, if necessary.
    // The allowable range is related to the PIVMIN.value parameter; see the
    // comments in DLARRD.  The preference for scaling small values
    // up is heuristic; we expect users' matrices not to be close to the
    // RMAX threshold.

    SCALE = ONE;
    TNRM = dlanst('M', N, D, E);
    if (TNRM > ZERO && TNRM < RMIN) {
      SCALE = RMIN / TNRM;
    } else if (TNRM > RMAX) {
      SCALE = RMAX / TNRM;
    }
    if (SCALE != ONE) {
      dscal(N, SCALE, D, 1);
      dscal(N - 1, SCALE, E, 1);
      TNRM = TNRM * SCALE;
      if (VALEIG) {
        // If eigenvalues in interval have to be found,
        // scale (WL.value, WU.value] accordingly
        WL.value = WL.value * SCALE;
        WU.value = WU.value * SCALE;
      }
    }

    // Compute the desired eigenvalues of the tridiagonal after splitting
    // into smaller subblocks if the corresponding off-diagonal elements
    // are small
    // THRESH is the splitting parameter for DLARRE
    // A negative THRESH forces the old splitting criterion based on the
    // size of the off-diagonal. A positive THRESH switches to splitting
    // which preserves relative accuracy.

    if (TRYRAC.value) {
      // Test whether the matrix warrants the more expensive relative approach.
      dlarrr(N, D, E, IINFO);
    } else {
      // The user does not care about relative accurately eigenvalues
      IINFO.value = -1;
    }
    // Set the splitting criterion
    if (IINFO.value == 0) {
      THRESH = EPS;
    } else {
      THRESH = -EPS;
      // relative accuracy is desired but T does not guarantee it
      TRYRAC.value = false;
    }

    if (TRYRAC.value) {
      // Copy original diagonal, needed to guarantee relative accuracy
      dcopy(N, D, 1, WORK(INDD), 1);
    }
    // Store the squares of the offdiagonal values of T
    for (J = 1; J <= N - 1; J++) {
      WORK[INDE2 + J - 1] = pow(E[J], 2).toDouble();
    }

    // Set the tolerance parameters for bisection
    if (!WANTZ) {
      // DLARRE computes the eigenvalues to full precision.
      RTOL1.value = FOUR * EPS;
      RTOL2.value = FOUR * EPS;
    } else {
      // DLARRE computes the eigenvalues to less than full precision.
      // DLARRV will refine the eigenvalue approximations, and we can
      // need less accurate initial bisection in DLARRE.
      // Note: these settings do only affect the subset case and DLARRE
      RTOL1.value = sqrt(EPS);
      RTOL2.value = max(sqrt(EPS) * 5.0e-3, FOUR * EPS);
    }
    dlarre(
        RANGE,
        N,
        WL,
        WU,
        IIL,
        IIU,
        D,
        E,
        WORK(INDE2),
        RTOL1.value,
        RTOL2.value,
        THRESH,
        NSPLIT,
        IWORK(IINSPL),
        M,
        W,
        WORK(INDERR),
        WORK(INDGP),
        IWORK(IINDBL),
        IWORK(IINDW),
        WORK(INDGRS),
        PIVMIN,
        WORK(INDWRK),
        IWORK(IINDWK),
        IINFO);
    if (IINFO.value != 0) {
      INFO.value = 10 + (IINFO.value).abs();
      return;
    }
    // Note that if RANGE != 'V', DLARRE computes bounds on the desired
    // part of the spectrum. All desired eigenvalues are contained in
    // (WL.value,WU.value]

    if (WANTZ) {
      // Compute the desired eigenvectors corresponding to the computed
      // eigenvalues

      dlarrv(
          N,
          WL.value,
          WU.value,
          D,
          E,
          PIVMIN.value,
          IWORK(IINSPL),
          M.value,
          1,
          M.value,
          MINRGP,
          RTOL1,
          RTOL2,
          W,
          WORK(INDERR),
          WORK(INDGP),
          IWORK(IINDBL),
          IWORK(IINDW),
          WORK(INDGRS),
          Z,
          LDZ,
          ISUPPZ,
          WORK(INDWRK),
          IWORK(IINDWK),
          IINFO);
      if (IINFO.value != 0) {
        INFO.value = 20 + (IINFO.value).abs();
        return;
      }
    } else {
      // DLARRE computes eigenvalues of the (shifted) root representation
      // DLARRV returns the eigenvalues of the unshifted matrix.
      // However, if the eigenvectors are not desired by the user, we need
      // to apply the corresponding shifts from DLARRE to obtain the
      // eigenvalues of the original matrix.
      for (J = 1; J <= M.value; J++) {
        ITMP.value = IWORK[IINDBL + J - 1];
        W[J] = W[J] + E[IWORK[IINSPL + ITMP.value - 1]];
      }
    }

    if (TRYRAC.value) {
      // Refine computed eigenvalues so that they are relatively accurate
      // with respect to the original matrix T.
      IBEGIN = 1;
      WBEGIN = 1;
      for (JBLK = 1; JBLK <= IWORK[IINDBL + M.value - 1]; JBLK++) {
        IEND = IWORK[IINSPL + JBLK - 1];
        IN = IEND - IBEGIN + 1;
        WEND = WBEGIN - 1;
        // check if any eigenvalues have to be refined in this block
        while (WEND < M.value) {
          if (IWORK[IINDBL + WEND] != JBLK) break;

          WEND++;
        }
        if (WEND < WBEGIN) {
          IBEGIN = IEND + 1;
          continue;
        }

        OFFSET = IWORK[IINDW + WBEGIN - 1] - 1;
        IFIRST = IWORK[IINDW + WBEGIN - 1];
        ILAST = IWORK[IINDW + WEND - 1];
        RTOL2.value = FOUR * EPS;
        dlarrj(
            IN,
            WORK(INDD + IBEGIN - 1),
            WORK(INDE2 + IBEGIN - 1),
            IFIRST,
            ILAST,
            RTOL2.value,
            OFFSET,
            W(WBEGIN),
            WORK(INDERR + WBEGIN - 1),
            WORK(INDWRK),
            IWORK(IINDWK),
            PIVMIN.value,
            TNRM,
            IINFO);
        IBEGIN = IEND + 1;
        WBEGIN = WEND + 1;
      }
    }

    // If matrix was scaled, then rescale eigenvalues appropriately.

    if (SCALE != ONE) {
      dscal(M.value, ONE / SCALE, W, 1);
    }
  }

  // If eigenvalues are not in increasing order, then sort them,
  // possibly along with eigenvectors.

  if (NSPLIT.value > 1 || N == 2) {
    if (!WANTZ) {
      dlasrt('I', M.value, W, IINFO);
      if (IINFO.value != 0) {
        INFO.value = 3;
        return;
      }
    } else {
      for (J = 1; J <= M.value - 1; J++) {
        I = 0;
        TMP = W[J];
        for (JJ = J + 1; JJ <= M.value; JJ++) {
          if (W[JJ] < TMP) {
            I = JJ;
            TMP = W[JJ];
          }
        }
        if (I != 0) {
          W[I] = W[J];
          W[J] = TMP;
          if (WANTZ) {
            dswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
            ITMP.value = ISUPPZ[2 * I - 1];
            ISUPPZ[2 * I - 1] = ISUPPZ[2 * J - 1];
            ISUPPZ[2 * J - 1] = ITMP.value;
            ITMP.value = ISUPPZ[2 * I];
            ISUPPZ[2 * I] = ISUPPZ[2 * J];
            ISUPPZ[2 * J] = ITMP.value;
          }
        }
      }
    }
  }

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
