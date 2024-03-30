import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dstev(
  final String JOBZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool WANTZ;
  int IMAX, ISCALE;
  double BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA = 0, SMLNUM, TNRM;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -6;
  }

  if (INFO.value != 0) {
    xerbla('DSTEV ', -INFO.value);
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
  // eigenvectors, call DSTEQR.

  if (!WANTZ) {
    dsterf(N, D, E, INFO);
  } else {
    dsteqr('I', N, D, E, Z, LDZ, WORK, INFO);
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

  if (ISCALE == 1) {
    if (INFO.value == 0) {
      IMAX = N;
    } else {
      IMAX = INFO.value - 1;
    }
    dscal(IMAX, ONE / SIGMA, D, 1);
  }
}
