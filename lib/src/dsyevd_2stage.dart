// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/dsytrd_2stage.dart';
import 'package:dart_lapack/src/ilaenv2stage.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsyevd_2stage(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> W_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final W = W_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int INDE,
      INDTAU,
      // INDWK2,
      INDWRK,
      ISCALE,
      LIWMIN = 0,
      LLWORK,
      // LLWRK2,
      LWMIN = 0,
      LHTRD = 0,
      LWTRD,
      KD,
      IB,
      INDHOUS;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');
  LQUERY = (LWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (!(lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(LOWER || lsame(UPLO, 'U'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }

  if (INFO.value == 0) {
    if (N <= 1) {
      LIWMIN = 1;
      LWMIN = 1;
    } else {
      KD = ilaenv2stage(1, 'DSYTRD_2STAGE', JOBZ, N, -1, -1, -1);
      IB = ilaenv2stage(2, 'DSYTRD_2STAGE', JOBZ, N, KD, -1, -1);
      LHTRD = ilaenv2stage(3, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1);
      LWTRD = ilaenv2stage(4, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1);
      if (WANTZ) {
        LIWMIN = 3 + 5 * N;
        LWMIN = 1 + 6 * N + 2 * pow(N, 2).toInt();
      } else {
        LIWMIN = 1;
        LWMIN = 2 * N + 1 + LHTRD + LWTRD;
      }
    }
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -10;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSYEVD_2STAGE', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = A[1][1];
    if (WANTZ) A[1][1] = ONE;
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

  ANRM = dlansy('M', UPLO, N, A, LDA, WORK);
  ISCALE = 0;
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) dlascl(UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO);

  // Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.

  INDE = 1;
  INDTAU = INDE + N;
  INDHOUS = INDTAU + N;
  INDWRK = INDHOUS + LHTRD;
  LLWORK = LWORK - INDWRK + 1;
  // INDWK2 = INDWRK + N * N;
  // LLWRK2 = LWORK - INDWK2 + 1;

  dsytrd_2stage(JOBZ, UPLO, N, A, LDA, W, WORK(INDE), WORK(INDTAU),
      WORK(INDHOUS), LHTRD, WORK(INDWRK), LLWORK, IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
  // tridiagonal matrix, then call DORMTR to multiply it by the
  // Householder transformations stored in A.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    // Not available in this release, and argument checking should not
    // let it getting here
    return;
    // dstedc('I', N, W, WORK(INDE), WORK(INDWRK), N, WORK(INDWK2), LLWRK2, IWORK,
    //     LIWORK, INFO);
    // dormtr('L', UPLO, 'N', N, N, A, LDA, WORK(INDTAU), WORK(INDWRK), N,
    //     WORK(INDWK2), LLWRK2, IINFO);
    // dlacpy('A', N, N, WORK(INDWRK), N, A, LDA);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) dscal(N, ONE / SIGMA, W, 1);

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
