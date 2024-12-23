// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhbtrd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanhb.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zstedc.dart';

void zhbevd(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int IMAX, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');
  LQUERY = (LWORK == -1 || LIWORK == -1 || LRWORK == -1);

  INFO.value = 0;
  if (N <= 1) {
    LWMIN = 1;
    LRWMIN = 1;
    LIWMIN = 1;
  } else {
    if (WANTZ) {
      LWMIN = 2 * pow(N, 2).toInt();
      LRWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
      LIWMIN = 3 + 5 * N;
    } else {
      LWMIN = N;
      LRWMIN = N;
      LIWMIN = 1;
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
    WORK[1] = LWMIN.toComplex();
    RWORK[1] = LRWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -11;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -13;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -15;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHBEVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = AB[1][1].real;
    if (WANTZ) Z[1][1] = Complex.one;
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

  ANRM = zlanhb('M', UPLO, N, KD, AB, LDAB, RWORK);
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
      zlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    } else {
      zlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    }
  }

  // Call ZHBTRD to reduce Hermitian band matrix to tridiagonal form.

  INDE = 1;
  INDWRK = INDE + N;
  INDWK2 = 1 + N * N;
  LLWK2 = LWORK - INDWK2 + 1;
  LLRWK = LRWORK - INDWRK + 1;
  zhbtrd(JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK(INDE), Z, LDZ, WORK, IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    zstedc('I', N, W, RWORK(INDE), WORK.asMatrix(N), N, WORK(INDWK2), LLWK2,
        RWORK(INDWRK), LLRWK, IWORK, LIWORK, INFO);
    zgemm('N', 'N', N, N, N, Complex.one, Z, LDZ, WORK.asMatrix(N), N,
        Complex.zero, WORK(INDWK2).asMatrix(N), N);
    zlacpy('A', N, N, WORK(INDWK2).asMatrix(N), N, Z, LDZ);
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

  WORK[1] = LWMIN.toComplex();
  RWORK[1] = LRWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
