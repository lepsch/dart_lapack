// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dstein.dart';
import 'package:lapack/src/dstemr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dstevr(
  final String JOBZ,
  final String RANGE,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<int> ISUPPZ_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final ISUPPZ = ISUPPZ_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool ALLEIG, INDEIG, TEST, LQUERY, VALEIG, WANTZ;
  String ORDER;
  int I,
      IEEEOK,
      IMAX,
      INDIBL,
      INDIFL,
      INDISP,
      INDIWO,
      ISCALE,
      ITMP1,
      J,
      JJ,
      LIWMIN,
      LWMIN;
  double BIGNUM,
      EPS,
      RMAX,
      RMIN,
      SAFMIN,
      SIGMA = 0,
      SMLNUM,
      TMP1,
      TNRM,
      VLL = 0,
      VUU = 0;
  final TRYRAC = Box(false);
  final NSPLIT = Box(0);

  // Test the input parameters.

  IEEEOK = ilaenv(10, 'DSTEVR', 'N', 1, 2, 3, 4);

  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

  LQUERY = ((LWORK == -1) || (LIWORK == -1));
  LWMIN = max(1, 20 * N);
  LIWMIN = max(1, 10 * N);

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) INFO.value = -7;
    } else if (INDEIG) {
      if (IL < 1 || IL > max(1, N)) {
        INFO.value = -8;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -9;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) {
      INFO.value = -14;
    }
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -17;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -19;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSTEVR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  if (N == 1) {
    if (ALLEIG || INDEIG) {
      M.value = 1;
      W[1] = D[1];
    } else {
      if (VL < D[1] && VU >= D[1]) {
        M.value = 1;
        W[1] = D[1];
      }
    }
    if (WANTZ) Z[1][1] = ONE;
    return;
  }

  // Get machine constants.

  SAFMIN = dlamch('Safe minimum');
  EPS = dlamch('Precision');
  SMLNUM = SAFMIN / EPS;
  BIGNUM = ONE / SMLNUM;
  RMIN = sqrt(SMLNUM);
  RMAX = min(sqrt(BIGNUM), ONE / sqrt(sqrt(SAFMIN)));

  // Scale matrix to allowable range, if necessary.

  ISCALE = 0;
  if (VALEIG) {
    VLL = VL;
    VUU = VU;
  }

  TNRM = dlanst('M', N, D, E);
  if (TNRM > ZERO && TNRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / TNRM;
  } else if (TNRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / TNRM;
  }
  if (ISCALE == 1) {
    dscal(N, SIGMA, D, 1);
    dscal(N - 1, SIGMA, E, 1);
    if (VALEIG) {
      VLL = VL * SIGMA;
      VUU = VU * SIGMA;
    }
  }

  // Initialize indices into workspaces.  Note: These indices are used only
  // if DSTERF or DSTEMR fail.

  // IWORK[INDIBL:INDIBL+M-1] corresponds to IBLOCK in DSTEBZ and
  // stores the block indices of each of the M<=N eigenvalues.
  INDIBL = 1;
  // IWORK[INDISP:INDISP+NSPLIT-1] corresponds to ISPLIT in DSTEBZ and
  // stores the starting and finishing indices of each block.
  INDISP = INDIBL + N;
  // IWORK[INDIFL:INDIFL+N-1] stores the indices of eigenvectors
  // that corresponding to eigenvectors that fail to converge in
  // DSTEIN.  This information is discarded; if any fail, the driver
  // returns INFO > 0.
  INDIFL = INDISP + N;
  // INDIWO is the offset of the remaining integer workspace.
  INDIWO = INDISP + N;

  // If all eigenvalues are desired, then
  // call DSTERF or DSTEMR.  If this fails for some eigenvalue, then
  // try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  while (true) {
    if ((ALLEIG || TEST) && IEEEOK == 1) {
      dcopy(N - 1, E, 1, WORK, 1);
      if (!WANTZ) {
        dcopy(N, D, 1, W, 1);
        dsterf(N, W, WORK, INFO);
      } else {
        dcopy(N, D, 1, WORK(N + 1), 1);
        if (ABSTOL <= TWO * N * EPS) {
          TRYRAC.value = true;
        } else {
          TRYRAC.value = false;
        }
        dstemr(
            JOBZ,
            'A',
            N,
            WORK(N + 1),
            WORK,
            VL,
            VU,
            IL,
            IU,
            M,
            W,
            Z,
            LDZ,
            N,
            ISUPPZ,
            TRYRAC,
            WORK(2 * N + 1),
            LWORK - 2 * N,
            IWORK,
            LIWORK,
            INFO);
      }
      if (INFO.value == 0) {
        M.value = N;
        break;
      }
      INFO.value = 0;
    }

    // Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, NSPLIT, W,
        IWORK(INDIBL), IWORK(INDISP), WORK, IWORK(INDIWO), INFO);

    if (WANTZ) {
      dstein(N, D, E, M.value, W, IWORK(INDIBL), IWORK(INDISP), Z, LDZ, WORK,
          IWORK(INDIWO), IWORK(INDIFL), INFO);
    }

    break;
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) {
    if (INFO.value == 0) {
      IMAX = M.value;
    } else {
      IMAX = INFO.value - 1;
    }
    dscal(IMAX, ONE / SIGMA, W, 1);
  }

  // If eigenvalues are not in order, then sort them, along with
  // eigenvectors.

  if (WANTZ) {
    for (J = 1; J <= M.value - 1; J++) {
      I = 0;
      TMP1 = W[J];
      for (JJ = J + 1; JJ <= M.value; JJ++) {
        if (W[JJ] < TMP1) {
          I = JJ;
          TMP1 = W[JJ];
        }
      }

      if (I != 0) {
        ITMP1 = IWORK[I];
        W[I] = W[J];
        IWORK[I] = IWORK[J];
        W[J] = TMP1;
        IWORK[J] = ITMP1;
        dswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
      }
    }
  }

  // Causes problems with tests 19 & 20:
  // IF (wantz && INDEIG ) Z[1][1] = Z[1][1] / 1.002 + .002

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
