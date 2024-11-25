// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhptrd.dart';
import 'package:dart_lapack/src/zlanhp.dart';
import 'package:dart_lapack/src/zstedc.dart';
import 'package:dart_lapack/src/zupmtr.dart';

void zhpevd(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
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
  final AP = AP_.having();
  final Z = Z_.having(ld: LDZ);
  final W = W_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, WANTZ;
  int IMAX,
      INDE,
      INDRWK,
      INDTAU,
      INDWRK,
      ISCALE,
      LIWMIN = 0,
      LLRWK,
      LLWRK,
      LRWMIN = 0,
      LWMIN = 0;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LQUERY = (LWORK == -1 || LRWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(lsame(UPLO, 'L') || lsame(UPLO, 'U'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -7;
  }

  if (INFO.value == 0) {
    if (N <= 1) {
      LWMIN = 1;
      LIWMIN = 1;
      LRWMIN = 1;
    } else {
      if (WANTZ) {
        LWMIN = 2 * N;
        LRWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
        LIWMIN = 3 + 5 * N;
      } else {
        LWMIN = N;
        LRWMIN = N;
        LIWMIN = 1;
      }
    }
    WORK[1] = LWMIN.toComplex();
    RWORK[1] = LRWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -9;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -11;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHPEVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = AP[1].real;
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

  ANRM = zlanhp('M', UPLO, N, AP, RWORK);
  ISCALE = 0;
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) {
    zdscal((N * (N + 1)) ~/ 2, SIGMA, AP, 1);
  }

  // Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.

  INDE = 1;
  INDTAU = 1;
  INDRWK = INDE + N;
  INDWRK = INDTAU + N;
  LLWRK = LWORK - INDWRK + 1;
  LLRWK = LRWORK - INDRWK + 1;
  zhptrd(UPLO, N, AP, W, RWORK(INDE), WORK(INDTAU), IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // ZUPGTR to generate the orthogonal matrix, then call ZSTEDC.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    zstedc('I', N, W, RWORK(INDE), Z, LDZ, WORK(INDWRK), LLWRK, RWORK(INDRWK),
        LLRWK, IWORK, LIWORK, INFO);
    zupmtr('L', UPLO, 'N', N, N, AP, WORK(INDTAU), Z, LDZ, WORK(INDWRK), IINFO);
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
