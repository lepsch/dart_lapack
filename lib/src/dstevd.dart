import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dstedc.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dstevd(
  final String JOBZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, WANTZ;
  int ISCALE, LIWMIN, LWMIN;
  double BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM, TNRM;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  LQUERY = (LWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  LIWMIN = 1;
  LWMIN = 1;
  if (N > 1 && WANTZ) {
    LWMIN = 1 + 4 * N + pow(N, 2).toInt();
    LIWMIN = 3 + 5 * N;
  }

  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -6;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -10;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSTEVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
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

  ISCALE = 0;
  TNRM = dlanst('M', N, D, E);
  if (TNRM > ZERO && TNRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / TNRM;
  } else if (TNRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / TNRM;
  }
  if (ISCALE == 1) {
    dscal(N, SIGMA, D, 1);
    dscal(N - 1, SIGMA, E(1), 1);
  }

  // For eigenvalues only, call DSTERF.  For eigenvalues and
  // eigenvectors, call DSTEDC.

  if (!WANTZ) {
    dsterf(N, D, E, INFO);
  } else {
    dstedc('I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) dscal(N, ONE / SIGMA, D, 1);

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
