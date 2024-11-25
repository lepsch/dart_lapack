// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhetrd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zstedc.dart';
import 'package:dart_lapack/src/zunmtr.dart';

void zheevd(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> W_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final W = W_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int IMAX,
      INDE,
      INDRWK,
      INDTAU,
      INDWK2,
      INDWRK,
      ISCALE,
      LIOPT = 0,
      LIWMIN,
      LLRWK,
      LLWORK,
      LLWRK2,
      LOPT = 0,
      LROPT = 0,
      LRWMIN,
      LWMIN;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');
  LQUERY = (LWORK == -1 || LRWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
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
      LWMIN = 1;
      LRWMIN = 1;
      LIWMIN = 1;
      LOPT = LWMIN;
      LROPT = LRWMIN;
      LIOPT = LIWMIN;
    } else {
      if (WANTZ) {
        LWMIN = 2 * N + N * N;
        LRWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
        LIWMIN = 3 + 5 * N;
      } else {
        LWMIN = N + 1;
        LRWMIN = N;
        LIWMIN = 1;
      }
      LOPT = max(LWMIN, N + N * ilaenv(1, 'ZHETRD', UPLO, N, -1, -1, -1));
      LROPT = LRWMIN;
      LIOPT = LIWMIN;
    }
    WORK[1] = LOPT.toComplex();
    RWORK[1] = LROPT.toDouble();
    IWORK[1] = LIOPT;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -10;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHEEVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = A[1][1].real;
    if (WANTZ) A[1][1] = Complex.one;
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

  ANRM = zlanhe('M', UPLO, N, A, LDA, RWORK);
  ISCALE = 0;
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) zlascl(UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO);

  // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

  INDE = 1;
  INDTAU = 1;
  INDWRK = INDTAU + N;
  INDRWK = INDE + N;
  INDWK2 = INDWRK + N * N;
  LLWORK = LWORK - INDWRK + 1;
  LLWRK2 = LWORK - INDWK2 + 1;
  LLRWK = LRWORK - INDRWK + 1;
  zhetrd(UPLO, N, A, LDA, W, RWORK(INDE), WORK(INDTAU), WORK(INDWRK), LLWORK,
      IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
  // tridiagonal matrix, then call ZUNMTR to multiply it to the
  // Householder transformations represented as Householder vectors in
  // A.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    zstedc('I', N, W, RWORK(INDE), WORK(INDWRK).asMatrix(), N, WORK(INDWK2),
        LLWRK2, RWORK(INDRWK), LLRWK, IWORK, LIWORK, INFO);
    zunmtr('L', UPLO, 'N', N, N, A, LDA, WORK(INDTAU), WORK(INDWRK).asMatrix(),
        N, WORK(INDWK2), LLWRK2, IINFO);
    zlacpy('A', N, N, WORK(INDWRK).asMatrix(), N, A, LDA);
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

  WORK[1] = LOPT.toComplex();
  RWORK[1] = LROPT.toDouble();
  IWORK[1] = LIOPT;
}
