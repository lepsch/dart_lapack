import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbtrd.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlanhb.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zstein.dart';
import 'package:lapack/src/zsteqr.dart';

void zhbevx(
  final String JOBZ,
  final String RANGE,
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> Q_,
  final int LDQ,
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
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final W = W_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ;
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
      INDWRK,
      ISCALE,
      ITMP1,
      J,
      JJ;
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
      VLL,
      VUU;
  Complex CTMP1;
  final IINFO = Box(0), NSPLIT = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  ALLEIG = lsame(RANGE, 'A');
  VALEIG = lsame(RANGE, 'V');
  INDEIG = lsame(RANGE, 'I');
  LOWER = lsame(UPLO, 'L');

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(ALLEIG || VALEIG || INDEIG)) {
    INFO.value = -2;
  } else if (!(LOWER || lsame(UPLO, 'U'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (KD < 0) {
    INFO.value = -5;
  } else if (LDAB < KD + 1) {
    INFO.value = -7;
  } else if (WANTZ && LDQ < max(1, N)) {
    INFO.value = -9;
  } else {
    if (VALEIG) {
      if (N > 0 && VU <= VL) INFO.value = -11;
    } else if (INDEIG) {
      if (IL < 1 || IL > max(1, N)) {
        INFO.value = -12;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -13;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N)) INFO.value = -18;
  }

  if (INFO.value != 0) {
    xerbla('ZHBEVX', -INFO.value);
    return;
  }

  // Quick return if possible

  M.value = 0;
  if (N == 0) return;

  if (N == 1) {
    M.value = 1;
    if (LOWER) {
      CTMP1 = AB[1][1];
    } else {
      CTMP1 = AB[KD + 1][1];
    }
    TMP1 = CTMP1.real;
    if (VALEIG) {
      if (!(VL < TMP1 && VU >= TMP1)) M.value = 0;
    }
    if (M.value == 1) {
      W[1] = CTMP1.real;
      if (WANTZ) Z[1][1] = Complex.one;
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
  } else {
    VLL = ZERO;
    VUU = ZERO;
  }
  ANRM = zlanhb('M', UPLO, N, KD, AB, LDAB, RWORK);
  if (ANRM > ZERO && ANRM < RMIN) {
    ISCALE = 1;
    SIGMA = RMIN / ANRM;
  } else if (ANRM > RMAX) {
    ISCALE = 1;
    SIGMA = RMAX / ANRM;
  }
  if (ISCALE == 1) {
    if (LOWER) {
      zlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    } else {
      zlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO);
    }
    if (ABSTOL > 0) ABSTLL = ABSTOL * SIGMA;
    if (VALEIG) {
      VLL = VL * SIGMA;
      VUU = VU * SIGMA;
    }
  }

  // Call ZHBTRD to reduce Hermitian band matrix to tridiagonal form.

  INDD = 1;
  INDE = INDD + N;
  INDRWK = INDE + N;
  INDWRK = 1;
  zhbtrd(JOBZ, UPLO, N, KD, AB, LDAB, RWORK(INDD), RWORK(INDE), Q, LDQ,
      WORK(INDWRK), IINFO);

  // If all eigenvalues are desired and ABSTOL is less than or equal
  // to zero, then call DSTERF or ZSTEQR.  If this fails for some
  // eigenvalue, then try DSTEBZ.

  TEST = false;
  if (INDEIG) {
    if (IL == 1 && IU == N) {
      TEST = true;
    }
  }
  var success = false;
  if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
    dcopy(N, RWORK(INDD), 1, W, 1);
    INDEE = INDRWK + 2 * N;
    if (!WANTZ) {
      dcopy(N - 1, RWORK(INDE), 1, RWORK(INDEE), 1);
      dsterf(N, W, RWORK(INDEE), INFO);
    } else {
      zlacpy('A', N, N, Q, LDQ, Z, LDZ);
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
      success = true;
    }
    INFO.value = 0;
  }

  if (!success) {
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

      for (J = 1; J <= M.value; J++) {
        zcopy(N, Z(1, J).asArray(), 1, WORK(1), 1);
        zgemv('N', N, N, Complex.one, Q, LDQ, WORK, 1, Complex.zero,
            Z(1, J).asArray(), 1);
      }
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
}
