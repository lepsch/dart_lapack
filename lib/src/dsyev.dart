import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dorgtr.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/dsytrd.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsyev(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> W_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final W = W_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWKOPT = 0, NB;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');
  LQUERY = (LWORK == -1);

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
    NB = ilaenv(1, 'DSYTRD', UPLO, N, -1, -1, -1);
    LWKOPT = max(1, (NB + 2) * N);
    WORK[1] = LWKOPT.toDouble();

    if (LWORK < max(1, 3 * N - 1) && !LQUERY) INFO.value = -8;
  }

  if (INFO.value != 0) {
    xerbla('DSYEV', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    return;
  }

  if (N == 1) {
    W[1] = A[1][1];
    WORK[1] = 2;
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

  // Call DSYTRD to reduce symmetric matrix to tridiagonal form.

  INDE = 1;
  INDTAU = INDE + N;
  INDWRK = INDTAU + N;
  LLWORK = LWORK - INDWRK + 1;
  dsytrd(UPLO, N, A, LDA, W, WORK(INDE), WORK(INDTAU), WORK(INDWRK), LLWORK,
      IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // DORGTR to generate the orthogonal matrix, then call DSTEQR.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    dorgtr(UPLO, N, A, LDA, WORK(INDTAU), WORK(INDWRK), LLWORK, IINFO);
    dsteqr(JOBZ, N, W, WORK(INDE), A, LDA, WORK(INDTAU), INFO);
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

  // Set WORK(1) to optimal workspace size.

  WORK[1] = LWKOPT.toDouble();
}
