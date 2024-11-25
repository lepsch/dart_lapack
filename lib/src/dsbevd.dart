// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlansb.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dsbtrd.dart';
import 'package:dart_lapack/src/dstedc.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsbevd(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLWRK2, LWMIN;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');
  LQUERY = (LWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (N <= 1) {
    LIWMIN = 1;
    LWMIN = 1;
  } else {
    if (WANTZ) {
      LIWMIN = 3 + 5 * N;
      LWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
    } else {
      LIWMIN = 1;
      LWMIN = 2 * N;
    }
  }
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(LOWER || lsame(UPLO, 'U'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (KD < 0) {
    INFO.value = -4;
  } else if (LDAB < KD + 1) {
    INFO.value = -6;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -9;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -11;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSBEVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = AB[1][1];
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

  ANRM = dlansb('M', UPLO, N, KD, AB, LDAB, WORK);
  ISCALE = 0;
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) {
    if (LOWER) {
      dlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    } else {
      dlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    }
  }

  // Call DSBTRD to reduce symmetric band matrix to tridiagonal form.

  INDE = 1;
  INDWRK = INDE + N;
  INDWK2 = INDWRK + N * N;
  LLWRK2 = LWORK - INDWK2 + 1;
  dsbtrd(
      JOBZ, UPLO, N, KD, AB, LDAB, W, WORK(INDE), Z, LDZ, WORK(INDWRK), IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    dstedc('I', N, W, WORK(INDE), WORK(INDWRK).asMatrix(N), N, WORK(INDWK2),
        LLWRK2, IWORK, LIWORK, INFO);
    dgemm('N', 'N', N, N, N, ONE, Z, LDZ, WORK(INDWRK).asMatrix(N), N, ZERO,
        WORK(INDWK2).asMatrix(N), N);
    dlacpy('A', N, N, WORK(INDWK2).asMatrix(N), N, Z, LDZ);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) dscal(N, ONE / SIGMA, W, 1);

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
