import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrd.dart';
import 'package:lapack/src/zlansy.dart';
import 'package:lapack/src/zstein.dart';
import 'package:lapack/src/zstemr.dart';
import 'package:lapack/src/zunmtr.dart';

void zheevr(
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<int> ISUPPZ_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final W = W_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final ISUPPZ = ISUPPZ_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ;
  String ORDER;
  int I,
      IEEEOK,
      IMAX,
      INDIBL,
      INDIFL,
      INDISP,
      INDIWO,
      INDRD,
      INDRDD,
      INDRE,
      INDREE,
      INDRWK,
      INDTAU,
      INDWK,
      INDWKN,
      ISCALE,
      ITMP1,
      J,
      JJ,
      LIWMIN,
      LLWORK,
      LLRWORK,
      LLWRKN,
      LRWMIN,
      LWKOPT = 0,
      LWMIN,
      NB;
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
  final TRYRAC = Box(false);

  // Test the input parameters.

  IEEEOK = ilaenv(10, 'ZHEEVR', 'N', 1, 2, 3, 4);

  LOWER = lsame(UPLO, 'L');
  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');

  LQUERY = ((LWORK == -1) || (LRWORK == -1) || (LIWORK == -1));

  if (N <= 1) {
    LWMIN = 1;
    LRWMIN = 1;
    LIWMIN = 1;
  } else {
    LWMIN = 2 * N;
    LRWMIN = 24 * N;
    LIWMIN = 10 * N;
  }

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
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
    NB = ilaenv(1, 'ZHETRD', UPLO, N, -1, -1, -1);
    NB = max(NB, ilaenv(1, 'ZUNMTR', UPLO, N, -1, -1, -1));
    LWKOPT = max((NB + 1) * N, LWMIN);
    WORK[1] = LWKOPT.toComplex();
    RWORK[1] = LRWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -18;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -20;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -22;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZHEEVR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) {
    WORK[1] = Complex.one;
    return;
  }

  if (N == 1) {
    WORK[1] = Complex.one;
    if (ALLEIG || INDEIG) {
      M.value = 1;
      W[1] = A[1][1].real;
    } else {
      if (VL < A[1][1].real && VU >= A[1][1].real) {
        M.value = 1;
        W[1] = A[1][1].real;
      }
    }
    if (WANTZ) {
      Z[1][1] = Complex.one;
      ISUPPZ[1] = 1;
      ISUPPZ[2] = 1;
    }
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
  ANRM = zlansy('M', UPLO, N, A, LDA, RWORK);
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
        zdscal(N - J + 1, SIGMA, A(J, J).asArray(), 1);
      }
    } else {
      for (J = 1; J <= N; J++) {
        zdscal(J, SIGMA, A(1, J).asArray(), 1);
      }
    }
    if (ABSTOL > 0) ABSTLL = ABSTOL * SIGMA;
    if (VALEIG) {
      VLL = VL * SIGMA;
      VUU = VU * SIGMA;
    }
  }

  // Initialize indices into workspaces.  Note: The IWORK indices are
  // used only if DSTERF or ZSTEMR fail.

  // WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
  // elementary reflectors used in ZHETRD.
  INDTAU = 1;
  // INDWK is the starting offset of the remaining complex workspace,
  // and LLWORK is the remaining complex workspace size.
  INDWK = INDTAU + N;
  LLWORK = LWORK - INDWK + 1;

  // RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
  // entries.
  INDRD = 1;
  // RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
  // tridiagonal matrix from ZHETRD.
  INDRE = INDRD + N;
  // RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
  // -written by ZSTEMR (the DSTERF path copies the diagonal to W).
  INDRDD = INDRE + N;
  // RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
  // -written while computing the eigenvalues in DSTERF and ZSTEMR.
  INDREE = INDRDD + N;
  // INDRWK is the starting offset of the left-over real workspace, and
  // LLRWORK is the remaining workspace size.
  INDRWK = INDREE + N;
  LLRWORK = LRWORK - INDRWK + 1;

  // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
  // stores the block indices of each of the M<=N eigenvalues.
  INDIBL = 1;
  // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
  // stores the starting and finishing indices of each block.
  INDISP = INDIBL + N;
  // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
  // that corresponding to eigenvectors that fail to converge in
  // DSTEIN.  This information is discarded; if any fail, the driver
  // returns INFO > 0.
  INDIFL = INDISP + N;
  // INDIWO is the offset of the remaining integer workspace.
  INDIWO = INDIFL + N;

  // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

  zhetrd(UPLO, N, A, LDA, RWORK(INDRD), RWORK(INDRE), WORK(INDTAU), WORK(INDWK),
      LLWORK, IINFO);

  // If all eigenvalues are desired
  // then call DSTERF or ZSTEMR and ZUNMTR.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  var success = false;
  if ((ALLEIG || TEST) && (IEEEOK == 1)) {
    if (!WANTZ) {
      dcopy(N, RWORK(INDRD), 1, W, 1);
      dcopy(N - 1, RWORK(INDRE), 1, RWORK(INDREE), 1);
      dsterf(N, W, RWORK(INDREE), INFO);
    } else {
      dcopy(N - 1, RWORK(INDRE), 1, RWORK(INDREE), 1);
      dcopy(N, RWORK(INDRD), 1, RWORK(INDRDD), 1);

      if (ABSTOL <= TWO * N * EPS) {
        TRYRAC.value = true;
      } else {
        TRYRAC.value = false;
      }
      zstemr(
          JOBZ,
          'A',
          N,
          RWORK(INDRDD),
          RWORK(INDREE),
          VL,
          VU,
          IL,
          IU,
          M,
          W,
          Z,
          LDZ,
          N,
          ISUPPZ,
          TRYRAC,
          RWORK(INDRWK),
          LLRWORK,
          IWORK,
          LIWORK,
          INFO);

      // Apply unitary matrix used in reduction to tridiagonal
      // form to eigenvectors returned by ZSTEMR.

      if (WANTZ && INFO.value == 0) {
        INDWKN = INDWK;
        LLWRKN = LWORK - INDWKN + 1;
        zunmtr('L', UPLO, 'N', N, M.value, A, LDA, WORK(INDTAU), Z, LDZ,
            WORK(INDWKN), LLWRKN, IINFO);
      }
    }

    if (INFO.value == 0) {
      M.value = N;
      success = true;
    }
    INFO.value = 0;
  }

  if (!success) {
    // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
    // Also call DSTEBZ and ZSTEIN if ZSTEMR fails.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    dstebz(
        RANGE,
        ORDER,
        N,
        VLL,
        VUU,
        IL,
        IU,
        ABSTLL,
        RWORK(INDRD),
        RWORK(INDRE),
        M,
        NSPLIT,
        W,
        IWORK(INDIBL),
        IWORK(INDISP),
        RWORK(INDRWK),
        IWORK(INDIWO),
        INFO);

    if (WANTZ) {
      zstein(
          N,
          RWORK(INDRD),
          RWORK(INDRE),
          M.value,
          W,
          IWORK(INDIBL),
          IWORK(INDISP),
          Z,
          LDZ,
          RWORK(INDRWK),
          IWORK(INDIWO),
          IWORK(INDIFL),
          INFO);

      // Apply unitary matrix used in reduction to tridiagonal
      // form to eigenvectors returned by ZSTEIN.

      INDWKN = INDWK;
      LLWRKN = LWORK - INDWKN + 1;
      zunmtr('L', UPLO, 'N', N, M.value, A, LDA, WORK(INDTAU), Z, LDZ,
          WORK(INDWKN), LLWRKN, IINFO);
    }

    // If matrix was scaled, then rescale eigenvalues appropriately.
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
        zswap(N, Z(1, I).asArray(), 1, Z(1, J).asArray(), 1);
      }
    }
  }

  // Set WORK(1) to optimal workspace size.

  WORK[1] = LWKOPT.toComplex();
  RWORK[1] = LRWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
