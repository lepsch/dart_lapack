import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrd_2stage.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zstedc.dart';
import 'package:lapack/src/zunmtr.dart';

void zheevd_2stage(
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

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER, LQUERY, WANTZ;
  int IMAX,
      INDE,
      INDRWK,
      INDTAU,
      INDWK2,
      INDWRK,
      ISCALE,
      LIWMIN = 0,
      LLRWK,
      LLWORK,
      LLWRK2,
      LRWMIN = 0,
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
  LQUERY = (LWORK == -1 || LRWORK == -1 || LIWORK == -1);

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
    if (N <= 1) {
      LWMIN = 1;
      LRWMIN = 1;
      LIWMIN = 1;
    } else {
      KD = ilaenv2stage(1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1);
      IB = ilaenv2stage(2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1);
      LHTRD = ilaenv2stage(3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
      LWTRD = ilaenv2stage(4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
      if (WANTZ) {
        LWMIN = 2 * N + N * N;
        LRWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
        LIWMIN = 3 + 5 * N;
      } else {
        LWMIN = N + 1 + LHTRD + LWTRD;
        LRWMIN = N;
        LIWMIN = 1;
      }
    }
    WORK[1] = LWMIN.toComplex();
    RWORK[1] = LRWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -10;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHEEVD_2STAGE', -INFO.value);
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

  // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

  INDE = 1;
  INDRWK = INDE + N;
  LLRWK = LRWORK - INDRWK + 1;
  INDTAU = 1;
  INDHOUS = INDTAU + N;
  INDWRK = INDHOUS + LHTRD;
  LLWORK = LWORK - INDWRK + 1;
  INDWK2 = INDWRK + N * N;
  LLWRK2 = LWORK - INDWK2 + 1;

  zhetrd_2stage(JOBZ, UPLO, N, A, LDA, W, RWORK(INDE), WORK(INDTAU),
      WORK(INDHOUS), LHTRD, WORK(INDWRK), LLWORK, IINFO);

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

  WORK[1] = LWMIN.toComplex();
  RWORK[1] = LRWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
