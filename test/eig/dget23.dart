// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeevx.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import 'dget22.dart';

void dget23(
  final bool COMP,
  final String BALANC,
  final int JTYPE,
  final double THRESH,
  final Array<int> ISEED_,
  final Nout NOUNIT,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> H_,
  final Array<double> WR_,
  final Array<double> WI_,
  final Array<double> WR1_,
  final Array<double> WI1_,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final Matrix<double> LRE_,
  final int LDLRE,
  final Array<double> RCONDV_,
  final Array<double> RCNDV1_,
  final Array<double> RCDVIN_,
  final Array<double> RCONDE_,
  final Array<double> RCNDE1_,
  final Array<double> RCDEIN_,
  final Array<double> SCALE_,
  final Array<double> SCALE1_,
  final Array<double> RESULT_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final WR = WR_.having();
  final WI = WI_.having();
  final WR1 = WR1_.having();
  final WI1 = WI1_.having();
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final LRE = LRE_.having(ld: LDLRE);
  final RCONDV = RCONDV_.having();
  final RCNDV1 = RCNDV1_.having();
  final RCDVIN = RCDVIN_.having();
  final RCONDE = RCONDE_.having();
  final RCNDE1 = RCNDE1_.having();
  final RCDEIN = RCDEIN_.having();
  final SCALE = SCALE_.having();
  final SCALE1 = SCALE1_.having();
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8;
  bool BALOK, NOBAL;
  String SENSE;
  int I, ISENS, ISENSM, J, JJ, KMIN;
  double EPS,
      SMLNUM,
      TNRM,
      TOL,
      TOLIN,
      ULP,
      ULPINV,
      V,
      VIMIN,
      VMAX,
      VMX,
      VRMIN,
      VRMX,
      VTST;
  final DUM = Array<double>(1), RES = Array<double>(2);
  const SENS = ['N', 'V'];
  final IINFO = Box(0),
      ILO = Box(0),
      IHI = Box(0),
      IHI1 = Box(0),
      ILO1 = Box(0);
  final ABNRM = Box(0.0), ABNRM1 = Box(0.0);

  // Check for errors

  NOBAL = lsame(BALANC, 'N');
  BALOK =
      NOBAL || lsame(BALANC, 'P') || lsame(BALANC, 'S') || lsame(BALANC, 'B');
  INFO.value = 0;
  if (!BALOK) {
    INFO.value = -2;
  } else if (THRESH < ZERO) {
    INFO.value = -4;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO = -6;
  } else if (N < 0) {
    INFO.value = -7;
  } else if (LDA < 1 || LDA < N) {
    INFO.value = -9;
  } else if (LDVL < 1 || LDVL < N) {
    INFO.value = -16;
  } else if (LDVR < 1 || LDVR < N) {
    INFO.value = -18;
  } else if (LDLRE < 1 || LDLRE < N) {
    INFO.value = -20;
  } else if (LWORK < 3 * N || (COMP && LWORK < 6 * N + N * N)) {
    INFO.value = -31;
  }

  if (INFO.value != 0) {
    xerbla('DGET23', -INFO.value);
    return;
  }

  // Quick return if nothing to do

  for (I = 1; I <= 11; I++) {
    RESULT[I] = -ONE;
  }

  if (N == 0) return;

  // More Important constants

  ULP = dlamch('Precision');
  SMLNUM = dlamch('S');
  ULPINV = ONE / ULP;

  // Compute eigenvalues and eigenvectors, and test them

  if (LWORK >= 6 * N + N * N) {
    SENSE = 'B';
    ISENSM = 2;
  } else {
    SENSE = 'E';
    ISENSM = 1;
  }
  dlacpy('F', N, N, A, LDA, H, LDA);
  dgeevx(BALANC, 'V', 'V', SENSE, N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO,
      IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, IINFO);
  if (IINFO.value != 0) {
    RESULT[1] = ULPINV;
    _printResults(NOUNIT, JTYPE, 'DGEEVX4', IINFO.value, N, BALANC, ISEED);
    INFO.value = (IINFO.value).abs();
    return;
  }

  // Do Test (1)

  dget22('N', 'N', 'N', N, A, LDA, VR, LDVR, WR, WI, WORK, RES);
  RESULT[1] = RES[1];

  // Do Test (2)

  dget22('T', 'N', 'T', N, A, LDA, VL, LDVL, WR, WI, WORK, RES);
  RESULT[2] = RES[1];

  // Do Test (3)

  for (J = 1; J <= N; J++) {
    TNRM = ONE;
    if (WI[J] == ZERO) {
      TNRM = dnrm2(N, VR(1, J).asArray(), 1);
    } else if (WI[J] > ZERO) {
      TNRM = dlapy2(
          dnrm2(N, VR(1, J).asArray(), 1), dnrm2(N, VR(1, J + 1).asArray(), 1));
    }
    RESULT[3] = max(RESULT[3], min(ULPINV, (TNRM - ONE).abs() / ULP));
    if (WI[J] > ZERO) {
      VMX = ZERO;
      VRMX = ZERO;
      for (JJ = 1; JJ <= N; JJ++) {
        VTST = dlapy2(VR[JJ][J], VR[JJ][J + 1]);
        if (VTST > VMX) VMX = VTST;
        if (VR[JJ][J + 1] == ZERO && VR[JJ][J].abs() > VRMX) {
          VRMX = VR[JJ][J].abs();
        }
      }
      if (VRMX / VMX < ONE - TWO * ULP) RESULT[3] = ULPINV;
    }
  }

  // Do Test (4)

  for (J = 1; J <= N; J++) {
    TNRM = ONE;
    if (WI[J] == ZERO) {
      TNRM = dnrm2(N, VL(1, J).asArray(), 1);
    } else if (WI[J] > ZERO) {
      TNRM = dlapy2(
          dnrm2(N, VL(1, J).asArray(), 1), dnrm2(N, VL(1, J + 1).asArray(), 1));
    }
    RESULT[4] = max(RESULT[4], min(ULPINV, (TNRM - ONE).abs() / ULP));
    if (WI[J] > ZERO) {
      VMX = ZERO;
      VRMX = ZERO;
      for (JJ = 1; JJ <= N; JJ++) {
        VTST = dlapy2(VL[JJ][J], VL[JJ][J + 1]);
        if (VTST > VMX) VMX = VTST;
        if (VL[JJ][J + 1] == ZERO && VL[JJ][J].abs() > VRMX) {
          VRMX = VL[JJ][J].abs();
        }
      }
      if (VRMX / VMX < ONE - TWO * ULP) RESULT[4] = ULPINV;
    }
  }

  // Test for all options of computing condition numbers

  for (ISENS = 1; ISENS <= ISENSM; ISENS++) {
    SENSE = SENS[ISENS - 1];

    // Compute eigenvalues only, and test them

    dlacpy('F', N, N, A, LDA, H, LDA);
    dgeevx(
        BALANC,
        'N',
        'N',
        SENSE,
        N,
        H,
        LDA,
        WR1,
        WI1,
        DUM.asMatrix(1),
        1,
        DUM.asMatrix(1),
        1,
        ILO1,
        IHI1,
        SCALE1,
        ABNRM1,
        RCNDE1,
        RCNDV1,
        WORK,
        LWORK,
        IWORK,
        IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      _printResults(NOUNIT, JTYPE, 'DGEEVX4', IINFO.value, N, BALANC, ISEED);
      INFO.value = (IINFO.value).abs();
      continue;
    }

    // Do Test (5)

    for (J = 1; J <= N; J++) {
      if (WR[J] != WR1[J] || WI[J] != WI1[J]) RESULT[5] = ULPINV;
    }

    // Do Test (8)

    if (!NOBAL) {
      for (J = 1; J <= N; J++) {
        if (SCALE[J] != SCALE1[J]) RESULT[8] = ULPINV;
      }
      if (ILO.value != ILO1.value) RESULT[8] = ULPINV;
      if (IHI.value != IHI1.value) RESULT[8] = ULPINV;
      if (ABNRM.value != ABNRM1.value) RESULT[8] = ULPINV;
    }

    // Do Test (9)

    if (ISENS == 2 && N > 1) {
      for (J = 1; J <= N; J++) {
        if (RCONDV[J] != RCNDV1[J]) RESULT[9] = ULPINV;
      }
    }

    // Compute eigenvalues and right eigenvectors, and test them

    dlacpy('F', N, N, A, LDA, H, LDA);
    dgeevx(
        BALANC,
        'N',
        'V',
        SENSE,
        N,
        H,
        LDA,
        WR1,
        WI1,
        DUM.asMatrix(1),
        1,
        LRE,
        LDLRE,
        ILO1,
        IHI1,
        SCALE1,
        ABNRM1,
        RCNDE1,
        RCNDV1,
        WORK,
        LWORK,
        IWORK,
        IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      _printResults(NOUNIT, JTYPE, 'DGEEVX4', IINFO.value, N, BALANC, ISEED);
      INFO.value = (IINFO.value).abs();
      continue;
    }

    // Do Test (5) again

    for (J = 1; J <= N; J++) {
      if (WR[J] != WR1[J] || WI[J] != WI1[J]) RESULT[5] = ULPINV;
    }

    // Do Test (6)

    for (J = 1; J <= N; J++) {
      for (JJ = 1; JJ <= N; JJ++) {
        if (VR[J][JJ] != LRE[J][JJ]) RESULT[6] = ULPINV;
      }
    }

    // Do Test (8) again

    if (!NOBAL) {
      for (J = 1; J <= N; J++) {
        if (SCALE[J] != SCALE1[J]) RESULT[8] = ULPINV;
      }
      if (ILO.value != ILO1.value) RESULT[8] = ULPINV;
      if (IHI.value != IHI1.value) RESULT[8] = ULPINV;
      if (ABNRM.value != ABNRM1.value) RESULT[8] = ULPINV;
    }

    // Do Test (9) again

    if (ISENS == 2 && N > 1) {
      for (J = 1; J <= N; J++) {
        if (RCONDV[J] != RCNDV1[J]) RESULT[9] = ULPINV;
      }
    }

    // Compute eigenvalues and left eigenvectors, and test them

    dlacpy('F', N, N, A, LDA, H, LDA);
    dgeevx(
        BALANC,
        'V',
        'N',
        SENSE,
        N,
        H,
        LDA,
        WR1,
        WI1,
        LRE,
        LDLRE,
        DUM.asMatrix(1),
        1,
        ILO1,
        IHI1,
        SCALE1,
        ABNRM1,
        RCNDE1,
        RCNDV1,
        WORK,
        LWORK,
        IWORK,
        IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      _printResults(NOUNIT, JTYPE, 'DGEEVX4', IINFO.value, N, BALANC, ISEED);
      INFO.value = (IINFO.value).abs();
      continue;
    }

    // Do Test (5) again

    for (J = 1; J <= N; J++) {
      if (WR[J] != WR1[J] || WI[J] != WI1[J]) RESULT[5] = ULPINV;
    }

    // Do Test (7)

    for (J = 1; J <= N; J++) {
      for (JJ = 1; JJ <= N; JJ++) {
        if (VL[J][JJ] != LRE[J][JJ]) RESULT[7] = ULPINV;
      }
    }

    // Do Test (8) again

    if (!NOBAL) {
      for (J = 1; J <= N; J++) {
        if (SCALE[J] != SCALE1[J]) RESULT[8] = ULPINV;
      }
      if (ILO.value != ILO1.value) RESULT[8] = ULPINV;
      if (IHI.value != IHI1.value) RESULT[8] = ULPINV;
      if (ABNRM.value != ABNRM1.value) RESULT[8] = ULPINV;
    }

    // Do Test (9) again

    if (ISENS == 2 && N > 1) {
      for (J = 1; J <= N; J++) {
        if (RCONDV[J] != RCNDV1[J]) RESULT[9] = ULPINV;
      }
    }
  }

  // If COMP, compare condition numbers to precomputed ones

  if (COMP) {
    dlacpy('F', N, N, A, LDA, H, LDA);
    dgeevx('N', 'V', 'V', 'B', N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI,
        SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      _printResults(NOUNIT, JTYPE, 'DGEEVX5', IINFO.value, N, BALANC, ISEED,
          comp: true);
      INFO.value = (IINFO.value).abs();
    } else {
      // Sort eigenvalues and condition numbers lexicographically
      // to compare with inputs

      for (I = 1; I <= N - 1; I++) {
        KMIN = I;
        VRMIN = WR[I];
        VIMIN = WI[I];
        for (J = I + 1; J <= N; J++) {
          if (WR[J] < VRMIN) {
            KMIN = J;
            VRMIN = WR[J];
            VIMIN = WI[J];
          }
        }
        WR[KMIN] = WR[I];
        WI[KMIN] = WI[I];
        WR[I] = VRMIN;
        WI[I] = VIMIN;
        VRMIN = RCONDE[KMIN];
        RCONDE[KMIN] = RCONDE[I];
        RCONDE[I] = VRMIN;
        VRMIN = RCONDV[KMIN];
        RCONDV[KMIN] = RCONDV[I];
        RCONDV[I] = VRMIN;
      }

      // Compare condition numbers for eigenvectors
      // taking their condition numbers into account

      RESULT[10] = ZERO;
      EPS = max(EPSIN, ULP);
      V = max(N * EPS * ABNRM.value, SMLNUM);
      if (ABNRM.value == ZERO) V = ONE;
      for (I = 1; I <= N; I++) {
        if (V > RCONDV[I] * RCONDE[I]) {
          TOL = RCONDV[I];
        } else {
          TOL = V / RCONDE[I];
        }
        if (V > RCDVIN[I] * RCDEIN[I]) {
          TOLIN = RCDVIN[I];
        } else {
          TOLIN = V / RCDEIN[I];
        }
        TOL = max(TOL, SMLNUM / EPS);
        TOLIN = max(TOLIN, SMLNUM / EPS);
        if (EPS * (RCDVIN[I] - TOLIN) > RCONDV[I] + TOL) {
          VMAX = ONE / EPS;
        } else if (RCDVIN[I] - TOLIN > RCONDV[I] + TOL) {
          VMAX = (RCDVIN[I] - TOLIN) / (RCONDV[I] + TOL);
        } else if (RCDVIN[I] + TOLIN < EPS * (RCONDV[I] - TOL)) {
          VMAX = ONE / EPS;
        } else if (RCDVIN[I] + TOLIN < RCONDV[I] - TOL) {
          VMAX = (RCONDV[I] - TOL) / (RCDVIN[I] + TOLIN);
        } else {
          VMAX = ONE;
        }
        RESULT[10] = max(RESULT[10], VMAX);
      }

      // Compare condition numbers for eigenvalues
      // taking their condition numbers into account

      RESULT[11] = ZERO;
      for (I = 1; I <= N; I++) {
        if (V > RCONDV[I]) {
          TOL = ONE;
        } else {
          TOL = V / RCONDV[I];
        }
        if (V > RCDVIN[I]) {
          TOLIN = ONE;
        } else {
          TOLIN = V / RCDVIN[I];
        }
        TOL = max(TOL, SMLNUM / EPS);
        TOLIN = max(TOLIN, SMLNUM / EPS);
        if (EPS * (RCDEIN[I] - TOLIN) > RCONDE[I] + TOL) {
          VMAX = ONE / EPS;
        } else if (RCDEIN[I] - TOLIN > RCONDE[I] + TOL) {
          VMAX = (RCDEIN[I] - TOLIN) / (RCONDE[I] + TOL);
        } else if (RCDEIN[I] + TOLIN < EPS * (RCONDE[I] - TOL)) {
          VMAX = ONE / EPS;
        } else if (RCDEIN[I] + TOLIN < RCONDE[I] - TOL) {
          VMAX = (RCONDE[I] - TOL) / (RCDEIN[I] + TOLIN);
        } else {
          VMAX = ONE;
        }
        RESULT[11] = max(RESULT[11], VMAX);
      }
    }
  }
}

void _printResults(
  final Nout nout,
  final int JTYPE,
  final String s,
  final int info,
  final int n,
  final String BALANC,
  final Array<int> iseed, {
  final bool comp = false,
}) {
  if (!comp && JTYPE != 22) {
    nout.println(
        ' DGET23: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${JTYPE.i6}, BALANC = $BALANC, ISEED=(${iseed.i5(4, ',')})');
  } else {
    nout.println(
        ' DGET23: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, INPUT EXAMPLE NUMBER = ${iseed[1].i4}');
  }
}
