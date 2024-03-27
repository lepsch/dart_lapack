import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeesx.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import 'common.dart';
import 'dort01.dart';
import 'dslect.dart';

void dget24(
  final bool COMP,
  final int JTYPE,
  final double THRESH,
  final Array<int> ISEED_,
  final Nout NOUNIT,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> H_,
  final Matrix<double> HT_,
  final Array<double> WR_,
  final Array<double> WI_,
  final Array<double> WRT_,
  final Array<double> WIT_,
  final Array<double> WRTMP_,
  final Array<double> WITMP_,
  final Matrix<double> VS_,
  final int LDVS,
  final Matrix<double> VS1_,
  final double RCDEIN,
  final double RCDVIN,
  final int NSLCT,
  final Array<int> ISLCT_,
  final Array<double> RESULT_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final HT = HT_.having(ld: LDA);
  final WR = WR_.having();
  final WI = WI_.having();
  final WRT = WRT_.having();
  final WIT = WIT_.having();
  final WRTMP = WRTMP_.having();
  final WITMP = WITMP_.having();
  final VS = VS_.having(ld: LDVS);
  final VS1 = VS1_.having(ld: LDVS);
  final ISLCT = ISLCT_.having();
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const EPSIN = 5.9605e-8;
  String SORT;
  int I, ISORT, ITMP, J, KMIN, KNTEIG, LIWORK, RSUB;
  double ANORM,
      EPS,
      SMLNUM,
      TMP,
      TOL,
      TOLIN,
      ULP,
      ULPINV,
      V,
      VIMIN,
      VRMIN,
      WNORM;
  final IPNT = Array<int>(20);
  final IINFO = Box(0), SDIM = Box(0), SDIM1 = Box(0);
  final RCONDE = Box(0.0),
      RCONDV = Box(0.0),
      RCNDE1 = Box(0.0),
      RCNDV1 = Box(0.0);

  // Check for errors

  INFO.value = 0;
  if (THRESH < ZERO) {
    INFO.value = -3;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO.value = -5;
  } else if (N < 0) {
    INFO.value = -6;
  } else if (LDA < 1 || LDA < N) {
    INFO.value = -8;
  } else if (LDVS < 1 || LDVS < N) {
    INFO.value = -18;
  } else if (LWORK < 3 * N) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('DGET24', -INFO.value);
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

  // Perform tests (1)-(13)
  main:
  while (true) {
    sslct.SELOPT = 0;
    LIWORK = N * N;
    for (ISORT = 0; ISORT <= 1; ISORT++) {
      if (ISORT == 0) {
        SORT = 'N';
        RSUB = 0;
      } else {
        SORT = 'S';
        RSUB = 6;
      }

      // Compute Schur form and Schur vectors, and test them

      dlacpy('F', N, N, A, LDA, H, LDA);
      dgeesx('V', SORT, dslect, 'N', N, H, LDA, SDIM, WR, WI, VS, LDVS, RCONDE,
          RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[1 + RSUB] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX1', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        return;
      }
      if (ISORT == 0) {
        dcopy(N, WR, 1, WRTMP, 1);
        dcopy(N, WI, 1, WITMP, 1);
      }

      // Do Test (1) or Test (7)

      RESULT[1 + RSUB] = ZERO;
      for (J = 1; J <= N - 2; J++) {
        for (I = J + 2; I <= N; I++) {
          if (H[I][J] != ZERO) RESULT[1 + RSUB] = ULPINV;
        }
      }
      for (I = 1; I <= N - 2; I++) {
        if (H[I + 1][I] != ZERO && H[I + 2][I + 1] != ZERO) {
          RESULT[1 + RSUB] = ULPINV;
        }
      }
      for (I = 1; I <= N - 1; I++) {
        if (H[I + 1][I] != ZERO) {
          if (H[I][I] != H[I + 1][I + 1] ||
              H[I][I + 1] == ZERO ||
              sign(ONE, H[I + 1][I]) == sign(ONE, H[I][I + 1])) {
            RESULT[1 + RSUB] = ULPINV;
          }
        }
      }

      // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)

      // Copy A to VS1, used as workspace

      dlacpy(' ', N, N, A, LDA, VS1, LDVS);

      // Compute Q*H and store in HT.

      dgemm('No transpose', 'No transpose', N, N, N, ONE, VS, LDVS, H, LDA,
          ZERO, HT, LDA);

      // Compute A - Q*H*Q'

      dgemm('No transpose', 'Transpose', N, N, N, -ONE, HT, LDA, VS, LDVS, ONE,
          VS1, LDVS);

      ANORM = max(dlange('1', N, N, A, LDA, WORK), SMLNUM);
      WNORM = dlange('1', N, N, VS1, LDVS, WORK);

      if (ANORM > WNORM) {
        RESULT[2 + RSUB] = (WNORM / ANORM) / (N * ULP);
      } else {
        if (ANORM < ONE) {
          RESULT[2 + RSUB] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
        } else {
          RESULT[2 + RSUB] = min(WNORM / ANORM, N) / (N * ULP);
        }
      }

      // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )

      dort01('Columns', N, N, VS, LDVS, WORK, LWORK, RESULT.box(3 + RSUB));

      // Do Test (4) or Test (10)

      RESULT[4 + RSUB] = ZERO;
      for (I = 1; I <= N; I++) {
        if (H[I][I] != WR[I]) RESULT[4 + RSUB] = ULPINV;
      }
      if (N > 1) {
        if (H[2][1] == ZERO && WI[1] != ZERO) RESULT[4 + RSUB] = ULPINV;
        if (H[N][N - 1] == ZERO && WI[N] != ZERO) RESULT[4 + RSUB] = ULPINV;
      }
      for (I = 1; I <= N - 1; I++) {
        if (H[I + 1][I] != ZERO) {
          TMP = sqrt(H[I + 1][I].abs()) * sqrt(H[I][I + 1].abs());
          RESULT[4 + RSUB] = max(
              RESULT[4 + RSUB], (WI[I] - TMP).abs() / max(ULP * TMP, SMLNUM));
          RESULT[4 + RSUB] = max(RESULT[4 + RSUB],
              (WI[I + 1] + TMP).abs() / max(ULP * TMP, SMLNUM));
        } else if (I > 1) {
          if (H[I + 1][I] == ZERO && H[I][I - 1] == ZERO && WI[I] != ZERO) {
            RESULT[4 + RSUB] = ULPINV;
          }
        }
      }

      // Do Test (5) or Test (11)

      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('N', SORT, dslect, 'N', N, HT, LDA, SDIM, WRT, WIT, VS, LDVS,
          RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[5 + RSUB] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX2', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
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
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[6 + RSUB] = ULPINV;
      }

      // Do Test (13)

      if (ISORT == 1) {
        RESULT[13] = ZERO;
        KNTEIG = 0;
        for (I = 1; I <= N; I++) {
          if (dslect(WR[I], WI[I]) || dslect(WR[I], -WI[I])) {
            KNTEIG++;
          }
          if (I < N) {
            if ((dslect(WR[I + 1], WI[I + 1]) ||
                    dslect(WR[I + 1], -WI[I + 1])) &&
                (!(dslect(WR[I], WI[I]) || dslect(WR[I], -WI[I]))) &&
                IINFO.value != N + 2) RESULT[13] = ULPINV;
          }
        }
        if (SDIM.value != KNTEIG) RESULT[13] = ULPINV;
      }
    }

    // If there is enough workspace, perform tests (14) and (15)
    // as well as (10) through (13)

    if (LWORK >= N + (N * N) ~/ 2) {
      // Compute both RCONDE.value and RCONDV.value with VS

      SORT = 'S';
      RESULT[14] = ZERO;
      RESULT[15] = ZERO;
      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('V', SORT, dslect, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
          RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[14] = ULPINV;
        RESULT[15] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX3', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
      }

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute both RCONDE.value and RCONDV.value without VS, and compare

      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('N', SORT, dslect, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
          RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[14] = ULPINV;
        RESULT[15] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX4', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
      }

      // Perform tests (14) and (15)

      if (RCNDE1.value != RCONDE.value) RESULT[14] = ULPINV;
      if (RCNDV1.value != RCONDV.value) RESULT[15] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDE.value with VS, and compare

      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('V', SORT, dslect, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
          RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[14] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX5', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
      }

      // Perform test (14)

      if (RCNDE1.value != RCONDE.value) RESULT[14] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDE.value without VS, and compare

      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('N', SORT, dslect, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
          RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[14] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX6', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
      }

      // Perform test (14)

      if (RCNDE1.value != RCONDE.value) RESULT[14] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDV.value with VS, and compare

      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('V', SORT, dslect, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
          RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[15] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX7', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
      }

      // Perform test (15)

      if (RCNDV1.value != RCONDV.value) RESULT[15] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[10] = ULPINV;
        for (J = 1; J <= N; J++) {
          if (H[I][J] != HT[I][J]) RESULT[11] = ULPINV;
          if (VS[I][J] != VS1[I][J]) RESULT[12] = ULPINV;
        }
      }
      if (SDIM.value != SDIM1.value) RESULT[13] = ULPINV;

      // Compute RCONDV.value without VS, and compare

      dlacpy('F', N, N, A, LDA, HT, LDA);
      dgeesx('N', SORT, dslect, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
          RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
      if (IINFO.value != 0 && IINFO.value != N + 2) {
        RESULT[15] = ULPINV;
        _printReport(NOUNIT, JTYPE, 'DGEESX8', IINFO.value, N, ISEED);
        INFO.value = IINFO.value.abs();
        break main;
      }

      // Perform test (15)

      if (RCNDV1.value != RCONDV.value) RESULT[15] = ULPINV;

      // Perform tests (10), (11), (12), and (13)

      for (I = 1; I <= N; I++) {
        if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[10] = ULPINV;
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
    // First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
    // the logical function dslect selects the eigenvalues specified
    // by NSLCT and ISLCT.

    sslct.SELDIM = N;
    sslct.SELOPT = 1;
    EPS = max(ULP, EPSIN);
    for (I = 1; I <= N; I++) {
      IPNT[I] = I;
      sslct.SELVAL[I] = false;
      sslct.SELWR[I] = WRTMP[I];
      sslct.SELWI[I] = WITMP[I];
    }
    for (I = 1; I <= N - 1; I++) {
      KMIN = I;
      VRMIN = WRTMP[I];
      VIMIN = WITMP[I];
      for (J = I + 1; J <= N; J++) {
        if (WRTMP[J] < VRMIN) {
          KMIN = J;
          VRMIN = WRTMP[J];
          VIMIN = WITMP[J];
        }
      }
      WRTMP[KMIN] = WRTMP[I];
      WITMP[KMIN] = WITMP[I];
      WRTMP[I] = VRMIN;
      WITMP[I] = VIMIN;
      ITMP = IPNT[I];
      IPNT[I] = IPNT[KMIN];
      IPNT[KMIN] = ITMP;
    }
    for (I = 1; I <= NSLCT; I++) {
      sslct.SELVAL[IPNT[ISLCT[I]]] = true;
    }

    // Compute condition numbers

    dlacpy('F', N, N, A, LDA, HT, LDA);
    dgeesx('N', 'S', dslect, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS,
        RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO);
    if (IINFO.value != 0 && IINFO.value != N + 2) {
      RESULT[16] = ULPINV;
      RESULT[17] = ULPINV;
      _printReport(NOUNIT, JTYPE, 'DGEESX9', IINFO.value, N, ISEED, comp: COMP);
      INFO.value = IINFO.value.abs();
      return;
    }

    // Compare condition number for average of selected eigenvalues
    // taking its condition number into account

    ANORM = dlange('1', N, N, A, LDA, WORK);
    V = max(N * EPS * ANORM, SMLNUM);
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

void _printReport(
  final Nout nout,
  final int jtype,
  final String s,
  final int info,
  final int n,
  final Array<int> iseed, {
  bool comp = false,
}) {
  if (!comp && jtype != 22) {
    //  NOUNIT.println( 9999, 'DGEESX1', info, N, JTYPE, ISEED);
    nout.println(
        ' DGET24: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
  } else {
    nout.println(
        ' DGET24: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, INPUT EXAMPLE NUMBER = ${iseed[1].i4}');
  }
}
