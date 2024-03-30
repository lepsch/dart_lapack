import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrd.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zsteqr.dart';
import 'package:lapack/src/zungtr.dart';

void zheev(
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
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final W = W_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
    NB = ilaenv(1, 'ZHETRD', UPLO, N, -1, -1, -1);
    LWKOPT = max(1, (NB + 1) * N);
    WORK[1] = LWKOPT.toComplex();

    if (LWORK < max(1, 2 * N - 1) && !LQUERY) INFO.value = -8;
  }

  if (INFO.value != 0) {
    xerbla('ZHEEV ', -INFO.value);
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

  // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

  INDE = 1;
  INDTAU = 1;
  INDWRK = INDTAU + N;
  LLWORK = LWORK - INDWRK + 1;
  zhetrd(UPLO, N, A, LDA, W, RWORK(INDE), WORK(INDTAU), WORK(INDWRK), LLWORK,
      IINFO);

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

  WORK[1] = LWKOPT.toComplex();
}
