// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlansp.dart';
import 'package:dart_lapack/src/dopgtr.dart';
import 'package:dart_lapack/src/dopmtr.dart';
import 'package:dart_lapack/src/dsptrd.dart';
import 'package:dart_lapack/src/dstebz.dart';
import 'package:dart_lapack/src/dstein.dart';
import 'package:dart_lapack/src/dsteqr.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dspevx(
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final Array<double> AP_,
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
  final AP = AP_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
  String ORDER;
  int I,
      IMAX,
      INDD,
      INDE,
      INDEE,
      INDISP,
      INDIWO,
      INDTAU,
      INDWRK,
      ISCALE,
      ITMP1,
      J,
      JJ;
  double ABSTLL,
      ANRM,
      BIGNUM,
      EPS,
      RMAX,
      RMIN,
      SAFMIN,
      SIGMA = 0,
      SMLNUM,
      TMP1,
      VLL,
      VUU;
  final IINFO = Box(0), NSPLIT = Box(0);

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
  } else if (!(lsame(UPLO, 'L') || lsame(UPLO, 'U'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
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
    xerbla('DSPEVX', -INFO.value);
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  if (N == 1) {
    if (ALLEIG || INDEIG) {
      M.value = 1;
      W[1] = AP[1];
    } else {
      if (VL < AP[1] && VU >= AP[1]) {
        M.value = 1;
        W[1] = AP[1];
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
  ABSTLL = ABSTOL;
  if (VALEIG) {
    VLL = VL;
    VUU = VU;
  } else {
    VLL = ZERO;
    VUU = ZERO;
  }
  ANRM = dlansp('M', UPLO, N, AP, WORK);
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) {
    dscal((N * (N + 1)) ~/ 2, SIGMA, AP, 1);
    if (ABSTOL > 0) ABSTLL = ABSTOL * SIGMA;
    if (VALEIG) {
      VLL = VL * SIGMA;
      VUU = VU * SIGMA;
    }
  }

  // Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.

  INDTAU = 1;
  INDE = INDTAU + N;
  INDD = INDE + N;
  INDWRK = INDD + N;
  dsptrd(UPLO, N, AP, WORK(INDD), WORK(INDE), WORK(INDTAU), IINFO);

  // If all eigenvalues are desired and ABSTOL is less than or equal
  // to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails
  // for some eigenvalue, then try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  while (true) {
    if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
      dcopy(N, WORK(INDD), 1, W, 1);
      INDEE = INDWRK + 2 * N;
      if (!WANTZ) {
        dcopy(N - 1, WORK(INDE), 1, WORK(INDEE), 1);
        dsterf(N, W, WORK(INDEE), INFO);
      } else {
        dopgtr(UPLO, N, AP, WORK(INDTAU), Z, LDZ, WORK(INDWRK), IINFO);
        dcopy(N - 1, WORK(INDE), 1, WORK(INDEE), 1);
        dsteqr(JOBZ, N, W, WORK(INDEE), Z, LDZ, WORK(INDWRK), INFO);
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
    INDISP = 1 + N;
    INDIWO = INDISP + N;
    dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK(INDD), WORK(INDE), M,
        NSPLIT, W, IWORK(1), IWORK(INDISP), WORK(INDWRK), IWORK(INDIWO), INFO);

    if (WANTZ) {
      dstein(N, WORK(INDD), WORK(INDE), M.value, W, IWORK(1), IWORK(INDISP), Z,
          LDZ, WORK(INDWRK), IWORK(INDIWO), IFAIL, INFO);

      // Apply orthogonal matrix used in reduction to tridiagonal
      // form to eigenvectors returned by DSTEIN.

      dopmtr('L', UPLO, 'N', N, M.value, AP, WORK(INDTAU), Z, LDZ, WORK(INDWRK),
          IINFO);
    }

    // If matrix was scaled, then rescale eigenvalues appropriately.

    break;
  }
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
