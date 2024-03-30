import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeevx.dart';
import 'package:lapack/src/zlacpy.dart';

import 'zget22.dart';

void zget23(
  final bool COMP,
  final int ISRT,
  final String BALANC,
  final int JTYPE,
  final double THRESH,
  final Array<int> ISEED_,
  final Nout NOUNIT,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final Array<Complex> W_,
  final Array<Complex> W1_,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Matrix<Complex> LRE_,
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
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final LRE = LRE_.having(ld: LDLRE);
  final W = W_.having();
  final W1 = W1_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RCONDV = RCONDV_.having();
  final RCNDV1 = RCNDV1_.having();
  final RCDVIN = RCDVIN_.having();
  final RCONDE = RCONDE_.having();
  final RCNDE1 = RCNDE1_.having();
  final RCDEIN = RCDEIN_.having();
  final SCALE = SCALE_.having();
  final SCALE1 = SCALE1_.having();
  final RESULT = RESULT_.having();
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
      VMAX,
      VMX,
      VRICMP,
      VRIMIN,
      VRMX,
      VTST;
  Complex CTMP;
  final RES = Array<double>(2);
  final CDUM = Array<Complex>(1);
  const SENS = ['N', 'V'];
  final IHI = Box(0),
      ILO = Box(0),
      IHI1 = Box(0),
      ILO1 = Box(0),
      IINFO = Box(0);
  final ABNRM = Box(0.0), ABNRM1 = Box(0.0);

  // Check for errors

  NOBAL = lsame(BALANC, 'N');
  BALOK =
      NOBAL || lsame(BALANC, 'P') || lsame(BALANC, 'S') || lsame(BALANC, 'B');
  INFO.value = 0;
  if (ISRT != 0 && ISRT != 1) {
    INFO.value = -2;
  } else if (!BALOK) {
    INFO.value = -3;
  } else if (THRESH < ZERO) {
    INFO.value = -5;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO = -7;
  } else if (N < 0) {
    INFO.value = -8;
  } else if (LDA < 1 || LDA < N) {
    INFO.value = -10;
  } else if (LDVL < 1 || LDVL < N) {
    INFO.value = -15;
  } else if (LDVR < 1 || LDVR < N) {
    INFO.value = -17;
  } else if (LDLRE < 1 || LDLRE < N) {
    INFO.value = -19;
  } else if (LWORK < 2 * N || (COMP && LWORK < 2 * N + N * N)) {
    INFO.value = -30;
  }

  if (INFO.value != 0) {
    xerbla('ZGET23', -INFO.value);
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

  if (LWORK >= 2 * N + N * N) {
    SENSE = 'B';
    ISENSM = 2;
  } else {
    SENSE = 'E';
    ISENSM = 1;
  }
  zlacpy('F', N, N, A, LDA, H, LDA);
  zgeevx(BALANC, 'V', 'V', SENSE, N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI,
      SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO);
  if (IINFO.value != 0) {
    RESULT[1] = ULPINV;
    if (JTYPE != 22) {
      _print9998(NOUNIT, 'ZGEEVX1', IINFO.value, N, JTYPE, BALANC, ISEED);
    } else {
      _print9999(NOUNIT, 'ZGEEVX1', IINFO.value, N, ISEED[1]);
    }
    INFO.value = (IINFO.value).abs();
    return;
  }

  // Do Test (1)

  zget22('N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES);
  RESULT[1] = RES[1];

  // Do Test (2)

  zget22('C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES);
  RESULT[2] = RES[1];

  // Do Test (3)

  for (J = 1; J <= N; J++) {
    TNRM = dznrm2(N, VR(1, J).asArray(), 1);
    RESULT[3] = max(RESULT[3], min(ULPINV, (TNRM - ONE).abs() / ULP));
    VMX = ZERO;
    VRMX = ZERO;
    for (JJ = 1; JJ <= N; JJ++) {
      VTST = VR[JJ][J].abs();
      if (VTST > VMX) VMX = VTST;
      if (VR[JJ][J].imaginary == ZERO && VR[JJ][J].real.abs() > VRMX) {
        VRMX = VR[JJ][J].real.abs();
      }
    }
    if (VRMX / VMX < ONE - TWO * ULP) RESULT[3] = ULPINV;
  }

  // Do Test (4)

  for (J = 1; J <= N; J++) {
    TNRM = dznrm2(N, VL(1, J).asArray(), 1);
    RESULT[4] = max(RESULT[4], min(ULPINV, (TNRM - ONE).abs() / ULP));
    VMX = ZERO;
    VRMX = ZERO;
    for (JJ = 1; JJ <= N; JJ++) {
      VTST = VL[JJ][J].abs();
      if (VTST > VMX) VMX = VTST;
      if (VL[JJ][J].imaginary == ZERO && VL[JJ][J].real.abs() > VRMX) {
        VRMX = VL[JJ][J].real.abs();
      }
    }
    if (VRMX / VMX < ONE - TWO * ULP) RESULT[4] = ULPINV;
  }

  // Test for all options of computing condition numbers

  for (ISENS = 1; ISENS <= ISENSM; ISENS++) {
    SENSE = SENS[ISENS - 1];

    // Compute eigenvalues only, and test them

    zlacpy('F', N, N, A, LDA, H, LDA);
    zgeevx(
        BALANC,
        'N',
        'N',
        SENSE,
        N,
        H,
        LDA,
        W1,
        CDUM.asMatrix(),
        1,
        CDUM.asMatrix(),
        1,
        ILO1,
        IHI1,
        SCALE1,
        ABNRM1,
        RCNDE1,
        RCNDV1,
        WORK,
        LWORK,
        RWORK,
        IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      if (JTYPE != 22) {
        _print9998(NOUNIT, 'ZGEEVX2', IINFO.value, N, JTYPE, BALANC, ISEED);
      } else {
        _print9999(NOUNIT, 'ZGEEVX2', IINFO.value, N, ISEED[1]);
      }
      INFO.value = (IINFO.value).abs();
      continue;
    }

    // Do Test (5)

    for (J = 1; J <= N; J++) {
      if (W[J] != W1[J]) RESULT[5] = ULPINV;
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

    zlacpy('F', N, N, A, LDA, H, LDA);
    zgeevx(
        BALANC,
        'N',
        'V',
        SENSE,
        N,
        H,
        LDA,
        W1,
        CDUM.asMatrix(),
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
        RWORK,
        IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      if (JTYPE != 22) {
        _print9998(NOUNIT, 'ZGEEVX3', IINFO.value, N, JTYPE, BALANC, ISEED);
      } else {
        _print9999(NOUNIT, 'ZGEEVX3', IINFO.value, N, ISEED[1]);
      }
      INFO.value = (IINFO.value).abs();
      continue;
    }

    // Do Test (5) again

    for (J = 1; J <= N; J++) {
      if (W[J] != W1[J]) RESULT[5] = ULPINV;
    }

    // Do Test (6)

    for (J = 1; J <= N; J++) {
      for (JJ = 1; JJ <= N; JJ++) {
        if (VR(J, JJ) != LRE(J, JJ)) RESULT[6] = ULPINV;
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

    zlacpy('F', N, N, A, LDA, H, LDA);
    zgeevx(
        BALANC,
        'V',
        'N',
        SENSE,
        N,
        H,
        LDA,
        W1,
        LRE,
        LDLRE,
        CDUM.asMatrix(),
        1,
        ILO1,
        IHI1,
        SCALE1,
        ABNRM1,
        RCNDE1,
        RCNDV1,
        WORK,
        LWORK,
        RWORK,
        IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      if (JTYPE != 22) {
        _print9998(NOUNIT, 'ZGEEVX4', IINFO.value, N, JTYPE, BALANC, ISEED);
      } else {
        _print9999(NOUNIT, 'ZGEEVX4', IINFO.value, N, ISEED[1]);
      }
      INFO.value = (IINFO.value).abs();
      continue;
    }

    // Do Test (5) again

    for (J = 1; J <= N; J++) {
      if (W[J] != W1[J]) RESULT[5] = ULPINV;
    }

    // Do Test (7)

    for (J = 1; J <= N; J++) {
      for (JJ = 1; JJ <= N; JJ++) {
        if (VL(J, JJ) != LRE(J, JJ)) RESULT[7] = ULPINV;
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
    zlacpy('F', N, N, A, LDA, H, LDA);
    zgeevx('N', 'V', 'V', 'B', N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI,
        SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = ULPINV;
      _print9999(NOUNIT, 'ZGEEVX5', IINFO.value, N, ISEED[1]);
      INFO.value = (IINFO.value).abs();
      return;
    }

    // Sort eigenvalues and condition numbers lexicographically
    // to compare with inputs

    for (I = 1; I <= N - 1; I++) {
      KMIN = I;
      if (ISRT == 0) {
        VRIMIN = W[I].real;
      } else {
        VRIMIN = W[I].imaginary;
      }
      for (J = I + 1; J <= N; J++) {
        if (ISRT == 0) {
          VRICMP = W[J].real;
        } else {
          VRICMP = W[J].imaginary;
        }
        if (VRICMP < VRIMIN) {
          KMIN = J;
          VRIMIN = VRICMP;
        }
      }
      CTMP = W[KMIN];
      W[KMIN] = W[I];
      W[I] = CTMP;
      VRIMIN = RCONDE[KMIN];
      RCONDE[KMIN] = RCONDE[I];
      RCONDE[I] = VRIMIN;
      VRIMIN = RCONDV[KMIN];
      RCONDV[KMIN] = RCONDV[I];
      RCONDV[I] = VRIMIN;
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

void _print9999(Nout nout, String s, int info, int n, int i) {
  nout.println(
      ' ZGET23: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, INPUT EXAMPLE NUMBER = ${i.i4}');
}

void _print9998(Nout nout, String s, int info, int n, int jtype, String balanc,
    Array<int> iseed) {
  nout.println(
      ' ZGET23: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, BALANC = $balanc, ISEED=(${iseed.i5(4, ',')})');
}
