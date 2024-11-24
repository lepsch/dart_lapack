// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/ilaenv2stage.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhetrd_2stage.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zsteqr.dart';
import 'package:dart_lapack/src/zungtr.dart';

void zheev_2stage(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> W_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final W = W_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int IMAX,
      INDE,
      INDTAU,
      INDWRK,
      ISCALE,
      LLWORK,
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
  LQUERY = (LWORK == -1);

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
    KD = ilaenv2stage(1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1);
    IB = ilaenv2stage(2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1);
    LHTRD = ilaenv2stage(3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
    LWTRD = ilaenv2stage(4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
    LWMIN = N + LHTRD + LWTRD;
    WORK[1] = LWMIN.toComplex();

    if (LWORK < LWMIN && !LQUERY) INFO.value = -8;
  }

  if (INFO.value != 0) {
    xerbla('ZHEEV_2STAGE', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    return;
  }

  if (N == 1) {
    W[1] = A[1][1].real;
    WORK[1] = Complex.one;
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

  // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

  INDE = 1;
  INDTAU = 1;
  INDHOUS = INDTAU + N;
  INDWRK = INDHOUS + LHTRD;
  LLWORK = LWORK - INDWRK + 1;

  zhetrd_2stage(JOBZ, UPLO, N, A, LDA, W, RWORK(INDE), WORK(INDTAU),
      WORK(INDHOUS), LHTRD, WORK(INDWRK), LLWORK, IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // ZUNGTR to generate the unitary matrix, then call ZSTEQR.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    zungtr(UPLO, N, A, LDA, WORK(INDTAU), WORK(INDWRK), LLWORK, IINFO);
    INDWRK = INDE + N;
    zsteqr(JOBZ, N, W, RWORK(INDE), A, LDA, RWORK(INDWRK), INFO);
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

  // Set WORK(1) to optimal complex workspace size.

  WORK[1] = LWMIN.toComplex();
}
