import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeesx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';

import 'common.dart';
import 'zslect.dart';
import 'zunt01.dart';

void zget24(
  final bool COMP,
  final int JTYPE,
  final double THRESH,
  final Array<int> ISEED_,
  final Nout NOUNIT,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final Matrix<Complex> HT_,
  final Array<Complex> W_,
  final Array<Complex> WT_,
  final Array<Complex> WTMP_,
  final Matrix<Complex> VS_,
  final int LDVS,
  final Matrix<Complex> VS1_,
  final double RCDEIN,
  final double RCDVIN,
  final int NSLCT,
  final Array<int> ISLCT_,
  final int ISRT,
  final Array<double> RESULT_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final HT = HT_.having(ld: LDA);
  final VS = VS_.having(ld: LDVS);
  final VS1 = VS1_.having(ld: LDVS);
  final W = W_.having();
  final WT = WT_.having();
  final WTMP = WTMP_.having();
  final ISLCT = ISLCT_.having();
  final RESULT = RESULT_.having(length: 17);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final BWORK = BWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const EPSIN = 5.9605e-8;
  String SORT;
  int I, ISORT, ITMP, J, KMIN, KNTEIG, RSUB;
  double ANORM, EPS, SMLNUM, TOL, TOLIN, ULP, ULPINV, V, VRICMP, VRIMIN, WNORM;
  Complex CTMP;
  final IPNT = Array<int>(20);
  final SDIM = Box(0), SDIM1 = Box(0), IINFO = Box(0);
  final RCONDE = Box(0.0),
      RCONDV = Box(0.0),
      RCNDE1 = Box(0.0),
      RCNDV1 = Box(0.0);

  // Check for errors

  INFO.value = 0;
  if (THRESH < ZERO) {
    INFO.value = -3;
    // } else if (NOUNIT <= 0) {
    //   INFO.value = -5;
  } else if (N < 0) {
    INFO.value = -6;
  } else if (LDA < 1 || LDA < N) {
    INFO.value = -8;
  } else if (LDVS < 1 || LDVS < N) {
    INFO.value = -15;
  } else if (LWORK < 2 * N) {
    INFO.value = -24;
  }

  if (INFO.value != 0) {
    xerbla('ZGET24', -INFO.value);
    return;
  }

  // Quick return if nothing to do

  for (I = 1; I <= 17; I++) {
    RESULT[I] = -ONE;
  }

  if (N == 0) return;

  // Important constants

  SMLNUM = dlamch('Safe minimum');
  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  performTests:
  while (true) {
    // Perform tests (1)-(13)

    sslct.SELOPT = 0;
    for (ISORT = 0; ISORT <= 1; ISORT++) {
      if (ISORT == 0) {
        SORT = 'N';
        RSUB = 0;
      } else {
        SORT = 'S';
        RSUB = 6;
      }

      // Compute Schur form and Schur vectors, and test them

      zlacpy('F', N, N, A, LDA, H, LDA);
      zgeesx('V', SORT, zslect, 'N', N, H, LDA, SDIM, W, VS, LDVS, RCONDE,
          RCONDV, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[1 + RSUB] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX1', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX1', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        return;
      }
      if (ISORT == 0) {
        zcopy(N, W, 1, WTMP, 1);
      }

      // Do Test (1) or Test (7)

      RESULT[1 + RSUB] = ZERO;
      for (J = 1; J <= N - 1; J++) {
        for (I = J + 1; I <= N; I++) {
          if (H[I][J] != Complex.zero) RESULT[1 + RSUB] = ULPINV;
        }
      }

      // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)

      // Copy A to VS1, used as workspace

      zlacpy(' ', N, N, A, LDA, VS1, LDVS);

      // Compute Q*H and store in HT.

      zgemm('No transpose', 'No transpose', N, N, N, Complex.one, VS, LDVS, H,
          LDA, Complex.zero, HT, LDA);

      // Compute A - Q*H*Q'

      zgemm('No transpose', 'Conjugate transpose', N, N, N, -Complex.one, HT,
          LDA, VS, LDVS, Complex.one, VS1, LDVS);

      ANORM = max(zlange('1', N, N, A, LDA, RWORK), SMLNUM);
      WNORM = zlange('1', N, N, VS1, LDVS, RWORK);

      if (ANORM > WNORM) {
        RESULT[2 + RSUB] = (WNORM / ANORM) / (N * ULP);
      } else {
        if (ANORM < ONE) {
          RESULT[2 + RSUB] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
        } else {
          RESULT[2 + RSUB] = min(WNORM / ANORM, N.toDouble()) / (N * ULP);
        }
      }

      // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )

      zunt01('Columns', N, N, VS, LDVS, WORK, LWORK, RWORK, RESULT(3 + RSUB));

      // Do Test (4) or Test (10)

      RESULT[4 + RSUB] = ZERO;
      for (I = 1; I <= N; I++) {
        if (H[I][I] != W[I]) RESULT[4 + RSUB] = ULPINV;
      }

      // Do Test (5) or Test (11)

      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('N', SORT, zslect, 'N', N, HT, LDA, SDIM, WT, VS, LDVS, RCONDE,
          RCONDV, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[5 + RSUB] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX2', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX2', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      RESULT[5 + RSUB] = ZERO;
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= N; I++) {
          if (H[I][J] != HT[I][J]) RESULT[5 + RSUB] = ULPINV;
        }
      }

      // Do Test (6) or Test (12)

      RESULT[6 + RSUB] = ZERO;
      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[6 + RSUB] = ULPINV;
      }

      // Do Test (13)

      if (ISORT == 1) {
        RESULT[13] = ZERO;
        KNTEIG = 0;
        for (I = 1; I <= N; I++) {
          if (zslect(W[I])) KNTEIG++;
          if (I < N) {
            if (zslect(W[I + 1]) && (!zslect(W[I]))) RESULT[13] = ULPINV;
          }
        }
        if (SDIM.value != KNTEIG) RESULT[13] = ULPINV;
      }
    }

    // If there is enough workspace, perform tests (14) and (15)
    // as well as (10) through (13)

    if (LWORK >= (N * (N + 1)) / 2) {
      // Compute both RCONDE.value and RCONDV.value with VS

      SORT = 'S';
      RESULT[14] = ZERO;
      RESULT[15] = ZERO;
      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('V', SORT, zslect, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE,
          RCONDV, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[14] = ULPINV;
        RESULT[15] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX3', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX3', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute both RCONDE.value and RCONDV.value without VS, and compare

      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('N', SORT, zslect, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1,
          RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[14] = ULPINV;
        RESULT[15] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX4', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX4', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      // Perform tests (14) and (15)

      if (RCNDE1.value != RCONDE.value) RESULT[14] = ULPINV;
      if (RCNDV1.value != RCONDV.value) RESULT[15] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDE.value with VS, and compare

      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('V', SORT, zslect, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1,
          RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[14] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX5', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX5', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      // Perform test (14)

      if (RCNDE1.value != RCONDE.value) RESULT[14] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDE.value without VS, and compare

      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('N', SORT, zslect, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1,
          RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[14] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX6', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX6', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      // Perform test (14)

      if (RCNDE1.value != RCONDE.value) RESULT[14] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDV.value with VS, and compare

      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('V', SORT, zslect, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1,
          RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[15] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX7', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX7', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      // Perform test (15)

      if (RCNDV1.value != RCONDV.value) RESULT[15] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDV.value without VS, and compare

      zlacpy('F', N, N, A, LDA, HT, LDA);
      zgeesx('N', SORT, zslect, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1,
          RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO);
      if (IINFO.value != 0) {
        RESULT[15] = ULPINV;
        if (JTYPE != 22) {
          _print9998(NOUNIT, 'ZGEESX8', IINFO.value, N, JTYPE, ISEED);
        } else {
          _print9999(NOUNIT, 'ZGEESX8', IINFO.value, N, ISEED[1]);
        }
        INFO.value = (IINFO.value).abs();
        break performTests;
      }

      // Perform test (15)

      if (RCNDV1.value != RCONDV.value) RESULT[15] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (W[I] != WT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;
    }
    break;
  }

  // If there are precomputed reciprocal condition numbers, compare
  // computed values with them.

  if (COMP) {
    // First set up sslct.SELOPT, sslct.SELDIM, sslct.SELVAL, sslct.SELWR and sslct.SELWI so that
    // the logical function zslect selects the eigenvalues specified
    // by NSLCT, ISLCT and ISRT.

    sslct.SELDIM = N;
    sslct.SELOPT = 1;
    EPS = max(ULP, EPSIN);
    for (I = 1; I <= N; I++) {
      IPNT[I] = I;
      sslct.SELVAL[I] = false;
      sslct.SELWR[I] = (WTMP[I]).toDouble();
      sslct.SELWI[I] = WTMP[I].imaginary;
    }
    for (I = 1; I <= N - 1; I++) {
      KMIN = I;
      if (ISRT == 0) {
        VRIMIN = (WTMP[I]).toDouble();
      } else {
        VRIMIN = WTMP[I].imaginary;
      }
      for (J = I + 1; J <= N; J++) {
        if (ISRT == 0) {
          VRICMP = (WTMP[J]).toDouble();
        } else {
          VRICMP = WTMP[J].imaginary;
        }
        if (VRICMP < VRIMIN) {
          KMIN = J;
          VRIMIN = VRICMP;
        }
      }
      CTMP = WTMP[KMIN];
      WTMP[KMIN] = WTMP[I];
      WTMP[I] = CTMP;
      ITMP = IPNT[I];
      IPNT[I] = IPNT[KMIN];
      IPNT[KMIN] = ITMP;
    }
    for (I = 1; I <= NSLCT; I++) {
      sslct.SELVAL[IPNT[ISLCT[I]]] = true;
    }

    // Compute condition numbers

    zlacpy('F', N, N, A, LDA, HT, LDA);
    zgeesx('N', 'S', zslect, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE,
        RCONDV, WORK, LWORK, RWORK, BWORK, IINFO);
    if (IINFO.value != 0) {
      RESULT[16] = ULPINV;
      RESULT[17] = ULPINV;
      _print9999(NOUNIT, 'ZGEESX9', IINFO.value, N, ISEED[1]);
      INFO.value = (IINFO.value).abs();
      return;
    }

    // Compare condition number for average of selected eigenvalues
    // taking its condition number into account

    ANORM = zlange('1', N, N, A, LDA, RWORK);
    V = max(N.toDouble() * EPS * ANORM, SMLNUM);
    if (ANORM == ZERO) V = ONE;
    if (V > RCONDV.value) {
      TOL = ONE;
    } else {
      TOL = V / RCONDV.value;
    }
    if (V > RCDVIN) {
      TOLIN = ONE;
    } else {
      TOLIN = V / RCDVIN;
    }
    TOL = max(TOL, SMLNUM / EPS);
    TOLIN = max(TOLIN, SMLNUM / EPS);
    if (EPS * (RCDEIN - TOLIN) > RCONDE.value + TOL) {
      RESULT[16] = ULPINV;
    } else if (RCDEIN - TOLIN > RCONDE.value + TOL) {
      RESULT[16] = (RCDEIN - TOLIN) / (RCONDE.value + TOL);
    } else if (RCDEIN + TOLIN < EPS * (RCONDE.value - TOL)) {
      RESULT[16] = ULPINV;
    } else if (RCDEIN + TOLIN < RCONDE.value - TOL) {
      RESULT[16] = (RCONDE.value - TOL) / (RCDEIN + TOLIN);
    } else {
      RESULT[16] = ONE;
    }

    // Compare condition numbers for right invariant subspace
    // taking its condition number into account

    if (V > RCONDV.value * RCONDE.value) {
      TOL = RCONDV.value;
    } else {
      TOL = V / RCONDE.value;
    }
    if (V > RCDVIN * RCDEIN) {
      TOLIN = RCDVIN;
    } else {
      TOLIN = V / RCDEIN;
    }
    TOL = max(TOL, SMLNUM / EPS);
    TOLIN = max(TOLIN, SMLNUM / EPS);
    if (EPS * (RCDVIN - TOLIN) > RCONDV.value + TOL) {
      RESULT[17] = ULPINV;
    } else if (RCDVIN - TOLIN > RCONDV.value + TOL) {
      RESULT[17] = (RCDVIN - TOLIN) / (RCONDV.value + TOL);
    } else if (RCDVIN + TOLIN < EPS * (RCONDV.value - TOL)) {
      RESULT[17] = ULPINV;
    } else if (RCDVIN + TOLIN < RCONDV.value - TOL) {
      RESULT[17] = (RCONDV.value - TOL) / (RCDVIN + TOLIN);
    } else {
      RESULT[17] = ONE;
    }
  }
}

void _print9999(Nout nout, String s, int info, int n, int i) {
  nout.println(
      ' ZGET24: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, INPUT EXAMPLE NUMBER = ${i.i4}');
}

void _print9998(
    Nout nout, String s, int info, int n, int jtype, Array<int> iseed) {
  nout.println(
      ' ZGET24: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
