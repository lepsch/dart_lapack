import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dorgtr.dart';
import 'package:lapack/src/dormtr.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dstein.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/dsytrd_2stage.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsyevx_2stage(
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W,
  final Matrix<double> Z,
  final int LDZ,
  final Array<double> WORK,
  final int LWORK,
  final Array<int> IWORK,
  final Array<int> IFAIL,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ;
  String ORDER;
  int I,
      IMAX,
      INDD,
      INDE,
      INDEE,
      INDIBL = 0,
      INDISP,
      INDIWO,
      INDTAU,
      INDWKN,
      INDWRK,
      ISCALE,
      ITMP1,
      J,
      JJ,
      LLWORK,
      LLWRKN,
      LWMIN = 0,
      LHTRD = 0,
      LWTRD,
      KD,
      IB,
      INDHOUS;
  double ABSTLL,
      ANRM,
      BIGNUM,
      EPS,
      RMAX,
      RMIN,
      SAFMIN,
      SIGMA = 0,
      SMLNUM,
      TMP1,
      VLL = 0,
      VUU = 0;
  final IINFO = Box(0), NSPLIT = Box(0);

  // Test the input parameters.

  LOWER = lsame(UPLO, 'L');
  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');
  LQUERY = (LWORK == -1);

  INFO.value = 0;
  if (!(lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -2;
  } else if (!(LOWER || lsame(UPLO, 'U'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) INFO.value = -8;
    } else if (INDEIG) {
      if (IL < 1 || IL > max(1, N)) {
        INFO.value = -9;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -10;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) {
      INFO.value = -15;
    }
  }

  if (INFO.value == 0) {
    if (N <= 1) {
      LWMIN = 1;
      WORK[1] = LWMIN.toDouble();
    } else {
      KD = ilaenv2stage(1, 'DSYTRD_2STAGE', JOBZ, N, -1, -1, -1);
      IB = ilaenv2stage(2, 'DSYTRD_2STAGE', JOBZ, N, KD, -1, -1);
      LHTRD = ilaenv2stage(3, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1);
      LWTRD = ilaenv2stage(4, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1);
      LWMIN = max(8 * N, 3 * N + LHTRD + LWTRD);
      WORK[1] = LWMIN.toDouble();
    }

    if (LWORK < LWMIN && !LQUERY) INFO.value = -17;
  }

  if (INFO.value != 0) {
    xerbla('DSYEVX_2STAGE', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) {
    return;
  }

  if (N == 1) {
    if (ALLEIG || INDEIG) {
      M.value = 1;
      W[1] = A[1][1];
    } else {
      if (VL < A[1][1] && VU >= A[1][1]) {
        M.value = 1;
        W[1] = A[1][1];
      }
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
  RMAX = min(sqrt(BIGNUM), ONE / sqrt(sqrt(SAFMIN)));

  // Scale matrix to allowable range, if necessary.

  ISCALE = 0;
  ABSTLL = ABSTOL;
  if (VALEIG) {
    VLL = VL;
    VUU = VU;
  }
  ANRM = dlansy('M.value', UPLO, N, A, LDA, WORK);
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) {
    if (LOWER) {
      for (J = 1; J <= N; J++) {
        dscal(N - J + 1, SIGMA, A(J, J).asArray(), 1);
      }
    } else {
      for (J = 1; J <= N; J++) {
        dscal(J, SIGMA, A(1, J).asArray(), 1);
      }
    }
    if (ABSTOL > 0) ABSTLL = ABSTOL * SIGMA;
    if (VALEIG) {
      VLL = VL * SIGMA;
      VUU = VU * SIGMA;
    }
  }

  // Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.

  INDTAU = 1;
  INDE = INDTAU + N;
  INDD = INDE + N;
  INDHOUS = INDD + N;
  INDWRK = INDHOUS + LHTRD;
  LLWORK = LWORK - INDWRK + 1;

  dsytrd_2stage(
    JOBZ,
    UPLO,
    N,
    A,
    LDA,
    WORK[INDD],
    WORK[INDE],
    WORK[INDTAU],
    WORK[INDHOUS],
    LHTRD,
    WORK[INDWRK],
    LLWORK,
    IINFO,
  );

  // If all eigenvalues are desired and ABSTOL is less than or equal to
  // zero, then call DSTERF or DORGTR and SSTEQR.  If this fails for
  // some eigenvalue, then try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  while (true) {
    if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
      dcopy(N, WORK(INDD), 1, W, 1);
      INDEE = INDWRK + 2 * N;
      if (!WANTZ) {
        dcopy(N - 1, WORK(INDE), 1, WORK(INDEE), 1);
        dsterf(N, W, WORK(INDEE), INFO);
      } else {
        dlacpy('A', N, N, A, LDA, Z, LDZ);
        dorgtr(UPLO, N, Z, LDZ, WORK(INDTAU), WORK(INDWRK), LLWORK, IINFO);
        dcopy(N - 1, WORK(INDE), 1, WORK(INDEE), 1);
        dsteqr(JOBZ, N, W, WORK(INDEE), Z, LDZ, WORK(INDWRK), INFO);
        if (INFO.value == 0) {
          for (I = 1; I <= N; I++) {
            IFAIL[I] = 0;
          }
        }
      }
      if (INFO.value == 0) {
        M.value = N;
        break;
      }
      INFO.value = 0;
    }

    // Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    INDIBL = 1;
    INDISP = INDIBL + N;
    INDIWO = INDISP + N;
    dstebz(
      RANGE,
      ORDER,
      N,
      VLL,
      VUU,
      IL,
      IU,
      ABSTLL,
      WORK(INDD),
      WORK(INDE),
      M,
      NSPLIT,
      W,
      IWORK(INDIBL),
      IWORK(INDISP),
      WORK(INDWRK),
      IWORK(INDIWO),
      INFO,
    );

    if (WANTZ) {
      dstein(
        N,
        WORK(INDD),
        WORK(INDE),
        M.value,
        W,
        IWORK(INDIBL),
        IWORK(INDISP),
        Z,
        LDZ,
        WORK(INDWRK),
        IWORK(INDIWO),
        IFAIL,
        INFO,
      );

      // Apply orthogonal matrix used in reduction to tridiagonal
      // form to eigenvectors returned by DSTEIN.

      INDWKN = INDE;
      LLWRKN = LWORK - INDWKN + 1;
      dormtr(
        'L',
        UPLO,
        'N',
        N,
        M.value,
        A,
        LDA,
        WORK(INDTAU),
        Z,
        LDZ,
        WORK(INDWKN),
        LLWRKN,
        IINFO,
      );
    }

    // If matrix was scaled, then rescale eigenvalues appropriately.

    break;
  }
  if (ISCALE == 1) {
    if (INFO.value == 0) {
      IMAX = M.value;
    } else {
      IMAX = INFO.value - 1;
    }
    dscal(IMAX, ONE / SIGMA, W, 1);
  }

  // If eigenvalues are not in order, then sort them, along with
  // eigenvectors.

  if (WANTZ) {
    for (J = 1; J <= M.value - 1; J++) {
      I = 0;
      TMP1 = W[J];
      for (JJ = J + 1; JJ <= M.value; JJ++) {
        if (W[JJ] < TMP1) {
          I = JJ;
          TMP1 = W[JJ];
        }
      }

      if (I != 0) {
        ITMP1 = IWORK[INDIBL + I - 1];
        W[I] = W[J];
        IWORK[INDIBL + I - 1] = IWORK[INDIBL + J - 1];
        W[J] = TMP1;
        IWORK[INDIBL + J - 1] = ITMP1;
        dswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
        if (INFO.value != 0) {
          ITMP1 = IFAIL[I];
          IFAIL[I] = IFAIL[J];
          IFAIL[J] = ITMP1;
        }
      }
    }
  }

  // Set WORK[1] to optimal workspace size.

  WORK[1] = LWMIN.toDouble();
}
