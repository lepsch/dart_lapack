import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrd_2stage.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zstein.dart';
import 'package:lapack/src/zsteqr.dart';
import 'package:lapack/src/zungtr.dart';
import 'package:lapack/src/zunmtr.dart';

void zheevx_2stage(
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
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  final W = W_.having();

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
      INDIWK,
      INDRWK,
      INDTAU,
      INDWRK,
      ISCALE,
      ITMP1,
      J,
      JJ,
      LLWORK,
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
      WORK[1] = LWMIN.toComplex();
    } else {
      KD = ilaenv2stage(1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1);
      IB = ilaenv2stage(2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1);
      LHTRD = ilaenv2stage(3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
      LWTRD = ilaenv2stage(4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1);
      LWMIN = N + LHTRD + LWTRD;
      WORK[1] = LWMIN.toComplex();
    }

    if (LWORK < LWMIN && !LQUERY) INFO.value = -17;
  }

  if (INFO.value != 0) {
    xerbla('ZHEEVX_2STAGE', -INFO.value);
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
      W[1] = A[1][1].real;
    } else if (VALEIG) {
      if (VL < A[1][1].real && VU >= A[1][1].real) {
        M.value = 1;
        W[1] = A[1][1].real;
      }
    }
    if (WANTZ) Z[1][1] = Complex.one;
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
  ANRM = zlanhe('M', UPLO, N, A, LDA, RWORK);
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

  // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

  INDD = 1;
  INDE = INDD + N;
  INDRWK = INDE + N;
  INDTAU = 1;
  INDHOUS = INDTAU + N;
  INDWRK = INDHOUS + LHTRD;
  LLWORK = LWORK - INDWRK + 1;

  zhetrd_2stage(JOBZ, UPLO, N, A, LDA, RWORK(INDD), RWORK(INDE), WORK(INDTAU),
      WORK(INDHOUS), LHTRD, WORK(INDWRK), LLWORK, IINFO);

  // If all eigenvalues are desired and ABSTOL is less than or equal to
  // zero, then call DSTERF or ZUNGTR and ZSTEQR.  If this fails for
  // some eigenvalue, then try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  var succeess = false;
  if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
    dcopy(N, RWORK(INDD), 1, W, 1);
    INDEE = INDRWK + 2 * N;
    if (!WANTZ) {
      dcopy(N - 1, RWORK(INDE), 1, RWORK(INDEE), 1);
      dsterf(N, W, RWORK(INDEE), INFO);
    } else {
      zlacpy('A', N, N, A, LDA, Z, LDZ);
      zungtr(UPLO, N, Z, LDZ, WORK(INDTAU), WORK(INDWRK), LLWORK, IINFO);
      dcopy(N - 1, RWORK(INDE), 1, RWORK(INDEE), 1);
      zsteqr(JOBZ, N, W, RWORK(INDEE), Z, LDZ, RWORK(INDRWK), INFO);
      if (INFO.value == 0) {
        for (I = 1; I <= N; I++) {
          IFAIL[I] = 0;
        }
      }
    }
    if (INFO.value == 0) {
      M.value = N;
      succeess = true;
    }
    INFO.value = 0;
  }

  if (!succeess) {
    // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.

    if (WANTZ) {
      ORDER = 'B';
    } else {
      ORDER = 'E';
    }
    INDIBL = 1;
    INDISP = INDIBL + N;
    INDIWK = INDISP + N;
    dstebz(
        RANGE,
        ORDER,
        N,
        VLL,
        VUU,
        IL,
        IU,
        ABSTLL,
        RWORK(INDD),
        RWORK(INDE),
        M,
        NSPLIT,
        W,
        IWORK(INDIBL),
        IWORK(INDISP),
        RWORK(INDRWK),
        IWORK(INDIWK),
        INFO);

    if (WANTZ) {
      zstein(N, RWORK(INDD), RWORK(INDE), M.value, W, IWORK(INDIBL),
          IWORK(INDISP), Z, LDZ, RWORK(INDRWK), IWORK(INDIWK), IFAIL, INFO);

      // Apply unitary matrix used in reduction to tridiagonal
      // form to eigenvectors returned by ZSTEIN.

      zunmtr('L', UPLO, 'N', N, M.value, A, LDA, WORK(INDTAU), Z, LDZ,
          WORK(INDWRK), LLWORK, IINFO);
    }
  }

  // If matrix was scaled, then rescale eigenvalues appropriately.

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
        if (INFO.value != 0) {
          ITMP1 = IFAIL[I];
          IFAIL[I] = IFAIL[J];
          IFAIL[J] = ITMP1;
        }
      }
    }
  }

  // Set WORK(1) to optimal complex workspace size.

  WORK[1] = LWMIN.toComplex();
}
