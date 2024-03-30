import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dopmtr.dart';
import 'package:lapack/src/dsptrd.dart';
import 'package:lapack/src/dstedc.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dspevd(
    final String JOBZ,
    final String UPLO,
    final int N,
    final Array<double> AP_,
    final Array<double> W_,
    final Matrix<double> Z_,
    final int LDZ,
    final Array<double> WORK_,
    final int LWORK,
    final Array<int> IWORK_,
    final int LIWORK,
    final Box<int> INFO) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, WANTZ;
  int INDE, INDTAU, INDWRK, ISCALE, LIWMIN = 0, LLWORK, LWMIN = 0;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LQUERY = (LWORK == -1 || LIWORK == -1);

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

  if (INFO.value == 0) {
    if (N <= 1) {
      LIWMIN = 1;
      LWMIN = 1;
    } else {
      if (WANTZ) {
        LIWMIN = 3 + 5 * N;
        LWMIN = 1 + 6 * N + pow(N, 2).toInt();
      } else {
        LIWMIN = 1;
        LWMIN = 2 * N;
      }
    }
    IWORK[1] = LIWMIN;
    WORK[1] = LWMIN.toDouble();

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -9;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -11;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSPEVD', -INFO.value);
    return;
  } else if (LQUERY) {
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
  // DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
  // tridiagonal matrix, then call DOPMTR to multiply it by the
  // Householder transformations represented in AP.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    INDWRK = INDTAU + N;
    LLWORK = LWORK - INDWRK + 1;
    dstedc('I', N, W, WORK(INDE), Z, LDZ, WORK(INDWRK), LLWORK, IWORK, LIWORK,
        INFO);
    dopmtr('L', UPLO, 'N', N, N, AP, WORK(INDTAU), Z, LDZ, WORK(INDWRK), IINFO);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) dscal(N, ONE / SIGMA, W, 1);

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
