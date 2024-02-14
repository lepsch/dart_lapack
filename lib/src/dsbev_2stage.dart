import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansb.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/dsytrd_sb2st.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsbev_2stage(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  final W = W_.dim();
  final Z = Z_.dim(LDZ);
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, WANTZ, LQUERY;
  int IMAX,
      INDE,
      INDWRK,
      ISCALE,
      LLWORK,
      LWMIN = 0,
      LHTRD = 0,
      LWTRD,
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
  } else if (KD < 0) {
    INFO.value = -4;
  } else if (LDAB < KD + 1) {
    INFO.value = -6;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -9;
  }

  if (INFO.value == 0) {
    if (N <= 1) {
      LWMIN = 1;
      WORK[1] = LWMIN.toDouble();
    } else {
      IB = ilaenv2stage(2, 'DSYTRD_SB2ST', JOBZ, N, KD, -1, -1);
      LHTRD = ilaenv2stage(3, 'DSYTRD_SB2ST', JOBZ, N, KD, IB, -1);
      LWTRD = ilaenv2stage(4, 'DSYTRD_SB2ST', JOBZ, N, KD, IB, -1);
      LWMIN = N + LHTRD + LWTRD;
      WORK[1] = LWMIN.toDouble();
    }

    if (LWORK < LWMIN && !LQUERY) INFO.value = -11;
  }

  if (INFO.value != 0) {
    xerbla('DSBEV_2STAGE ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    if (LOWER) {
      W[1] = AB[1][1];
    } else {
      W[1] = AB[KD + 1][1];
    }
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

  ANRM = dlansb('M', UPLO, N, KD, AB, LDAB, WORK);
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
      dlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    } else {
      dlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    }
  }

  // Call DSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form.

  INDE = 1;
  INDHOUS = INDE + N;
  INDWRK = INDHOUS + LHTRD;
  LLWORK = LWORK - INDWRK + 1;

  dsytrd_sb2st('N', JOBZ, UPLO, N, KD, AB, LDAB, W, WORK(INDE), WORK(INDHOUS),
      LHTRD, WORK(INDWRK), LLWORK, IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    dsteqr(JOBZ, N, W, WORK(INDE), Z, LDZ, WORK(INDWRK), INFO);
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

  WORK[1] = LWMIN.toDouble();
}
