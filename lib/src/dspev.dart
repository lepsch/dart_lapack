// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlansp.dart';
import 'package:dart_lapack/src/dopgtr.dart';
import 'package:dart_lapack/src/dsptrd.dart';
import 'package:dart_lapack/src/dsteqr.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dspev(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool WANTZ;
  int IMAX, INDE, INDTAU, INDWRK, ISCALE;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(lsame(UPLO, 'U') || lsame(UPLO, 'L'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    xerbla('DSPEV', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = AP[1];
    if (WANTZ) Z[1][1] = ONE;
    return;
  }

  // Get machine constants.

  SAFMIN = dlamch('Safe minimum');
  EPS = dlamch('Precision');
  SMLNUM = SAFMIN / EPS;
  BIGNUM = ONE / SMLNUM;
  RMIN = sqrt(SMLNUM);
  RMAX = sqrt(BIGNUM);

  // Scale matrix to allowable range, if necessary.

  ANRM = dlansp('M', UPLO, N, AP, WORK);
  ISCALE = 0;
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) {
    dscal((N * (N + 1)) ~/ 2, SIGMA, AP, 1);
  }

  // Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.

  INDE = 1;
  INDTAU = INDE + N;
  dsptrd(UPLO, N, AP, W, WORK(INDE), WORK(INDTAU), IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // DOPGTR to generate the orthogonal matrix, then call DSTEQR.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    INDWRK = INDTAU + N;
    dopgtr(UPLO, N, AP, WORK(INDTAU), Z, LDZ, WORK(INDWRK), IINFO);
    dsteqr(JOBZ, N, W, WORK(INDE), Z, LDZ, WORK(INDTAU), INFO);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) {
    if (INFO.value == 0) {
      IMAX = N;
    } else {
      IMAX = INFO.value - 1;
    }
    dscal(IMAX, ONE / SIGMA, W, 1);
  }
}
