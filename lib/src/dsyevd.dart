import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dormtr.dart';
import 'package:lapack/src/dstedc.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/dsytrd.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

// >
// =====================================================================
void dsyevd(
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
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final W = W_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();
  const ZERO = 0.0, ONE = 1.0;

  bool LOWER, LQUERY, WANTZ;
  int INDE,
      INDTAU,
      INDWK2,
      INDWRK,
      ISCALE,
      LIOPT = 0,
      LIWMIN,
      LLWORK,
      LLWRK2,
      LOPT = 0,
      LWMIN;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');
  LQUERY = (LWORK == -1 || LIWORK == -1);

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
      LIWMIN = 1;
      LWMIN = 1;
      LOPT = LWMIN;
      LIOPT = LIWMIN;
    } else {
      if (WANTZ) {
        LIWMIN = 3 + 5 * N;
        LWMIN = 1 + 6 * N + 2 * pow(N, 2).toInt();
      } else {
        LIWMIN = 1;
        LWMIN = 2 * N + 1;
      }
      LOPT = max(LWMIN, 2 * N + N * ilaenv(1, 'DSYTRD', UPLO, N, -1, -1, -1));
      LIOPT = LIWMIN;
    }
    WORK[1] = LOPT.toDouble();
    IWORK[1] = LIOPT;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -10;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSYEVD', -INFO.value);
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

  // Call DSYTRD to reduce symmetric matrix to tridiagonal form.

  INDE = 1;
  INDTAU = INDE + N;
  INDWRK = INDTAU + N;
  LLWORK = LWORK - INDWRK + 1;
  INDWK2 = INDWRK + N * N;
  LLWRK2 = LWORK - INDWK2 + 1;

  dsytrd(UPLO, N, A, LDA, W, WORK(INDE), WORK(INDTAU), WORK(INDWRK), LLWORK,
      IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
  // tridiagonal matrix, then call DORMTR to multiply it by the
  // Householder transformations stored in A.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    dstedc('I', N, W, WORK(INDE), WORK(INDWRK).asMatrix(N), N, WORK(INDWK2),
        LLWRK2, IWORK, LIWORK, INFO);
    dormtr('L', UPLO, 'N', N, N, A, LDA, WORK(INDTAU), WORK(INDWRK).asMatrix(N),
        N, WORK(INDWK2), LLWRK2, IINFO);
    dlacpy('A', N, N, WORK(INDWRK).asMatrix(N), N, A, LDA);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) dscal(N, ONE / SIGMA, W, 1);

  WORK[1] = LOPT.toDouble();
  IWORK[1] = LIOPT;
}
