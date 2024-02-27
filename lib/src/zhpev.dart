import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhptrd.dart';
import 'package:lapack/src/zlanhp.dart';
import 'package:lapack/src/zsteqr.dart';
import 'package:lapack/src/zupgtr.dart';

void zhpev(
  final String JOBZ,
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.dim(LDZ);
  final AP = AP_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final W = W_.dim();
  const ZERO = 0.0, ONE = 1.0;
  bool WANTZ;
  int IMAX, INDE, INDRWK, INDTAU, INDWRK, ISCALE;
  double ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');

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

  if (INFO.value != 0) {
    xerbla('ZHPEV ', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    W[1] = AP[1].toDouble();
    RWORK[1] = 1;
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
  zhptrd(UPLO, N, AP, W, RWORK(INDE), WORK(INDTAU), IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, first call
  // ZUPGTR to generate the orthogonal matrix, then call ZSTEQR.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    INDWRK = INDTAU + N;
    zupgtr(UPLO, N, AP, WORK(INDTAU), Z, LDZ, WORK(INDWRK), IINFO);
    INDRWK = INDE + N;
    zsteqr(JOBZ, N, W, RWORK(INDE), Z, LDZ, RWORK(INDRWK), INFO);
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
