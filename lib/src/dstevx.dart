// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlanst.dart';
import 'package:dart_lapack/src/dstebz.dart';
import 'package:dart_lapack/src/dstein.dart';
import 'package:dart_lapack/src/dsteqr.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dstevx(
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
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
  String ORDER;
  int I, IMAX, INDISP, INDIWO, INDWRK, ISCALE, ITMP1, J, JJ;
  double BIGNUM,
      EPS,
      RMAX,
      RMIN,
      SAFMIN,
      SIGMA = 0,
      SMLNUM,
      TMP1,
      TNRM,
      VLL,
      VUU;
  final NSPLIT = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

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
    if (LDZ < 1 || (WANTZ && LDZ < N)) INFO.value = -14;
  }

  if (INFO.value != 0) {
    xerbla('DSTEVX', -INFO.value);
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
  } else {
    VLL = ZERO;
    VUU = ZERO;
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

  // If all eigenvalues are desired and ABSTOL is less than zero, then
  // call DSTERF or SSTEQR.  If this fails for some eigenvalue, then
  // try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  while (true) {
    if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
      dcopy(N, D, 1, W, 1);
      dcopy(N - 1, E, 1, WORK, 1);
      INDWRK = N + 1;
      if (!WANTZ) {
        dsterf(N, W, WORK, INFO);
      } else {
        dsteqr('I', N, W, WORK, Z, LDZ, WORK(INDWRK), INFO);
        if (INFO.value == 0) {
          for (I = 1; I <= N; I++) {
            IFAIL[I] = 0;
          }
        }
      }
      if (INFO.value == 0) {
        M.value = N;
        break;
      }
      INFO.value = 0;
    }

    // Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    INDWRK = 1;
    INDISP = 1 + N;
    INDIWO = INDISP + N;
    dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IWORK,
        IWORK(INDISP), WORK(INDWRK), IWORK(INDIWO), INFO);

    if (WANTZ) {
      dstein(N, D, E, M.value, W, IWORK(1), IWORK(INDISP), Z, LDZ, WORK(INDWRK),
          IWORK(INDIWO), IFAIL, INFO);
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
        ITMP1 = IWORK[1 + I - 1];
        W[I] = W[J];
        IWORK[1 + I - 1] = IWORK[1 + J - 1];
        W[J] = TMP1;
        IWORK[1 + J - 1] = ITMP1;
        dswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
        if (INFO.value != 0) {
          ITMP1 = IFAIL[I];
          IFAIL[I] = IFAIL[J];
          IFAIL[J] = ITMP1;
        }
      }
    }
  }
}
