import 'dart:math';

import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgehrd.dart';
import 'package:lapack/src/zhseqr.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/ztrsen.dart';
import 'package:lapack/src/zunghr.dart';

import 'zhst01.dart';

Future<void> zget38(
  final Array<double> RMAX_,
  final Array<int> LMAX_,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final Nin NIN,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RMAX = RMAX_.dim(3);
  final LMAX = LMAX_.dim(3);
  final NINFO = NINFO_.dim(3);

  const LDT = 20, LWORK = 2 * LDT * (10 + LDT);
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8;
  int I, ISCL, ISRT, ITMP, J, KMIN, N, NDIM;
  double BIGNUM, EPS, SEPIN, SIN, SMLNUM, TNRM, TOL, TOLIN, V, VMAX, VMIN, VMUL;
  final SELECT = Array<bool>(LDT);
  final IPNT = Array<int>(LDT), ISELEC = Array<int>(LDT);
  final RESULT = Array<double>(2),
      RWORK = Array<double>(LDT),
      VAL = Array<double>(3),
      WSRT = Array<double>(LDT);
  final W = Array<Complex>(LDT),
      WORK = Array<Complex>(LWORK),
      WTMP = Array<Complex>(LDT);
  final Q = Matrix<Complex>(LDT, LDT),
      QSAV = Matrix<Complex>(LDT, LDT),
      QTMP = Matrix<Complex>(LDT, LDT),
      T = Matrix<Complex>(LDT, LDT),
      TMP = Matrix<Complex>(LDT, LDT),
      TSAV = Matrix<Complex>(LDT, LDT),
      TSAV1 = Matrix<Complex>(LDT, LDT),
      TTMP = Matrix<Complex>(LDT, LDT);
  final INFO = Box(0), M = Box(0);
  final S = Box(0.0), SEP = Box(0.0), STMP = Box(0.0), SEPTMP = Box(0.0);

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // EPSIN = 2**(-24) = precision to which input data computed

  EPS = max(EPS, EPSIN);
  RMAX[1] = ZERO;
  RMAX[2] = ZERO;
  RMAX[3] = ZERO;
  LMAX[1] = 0;
  LMAX[2] = 0;
  LMAX[3] = 0;
  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;
  NINFO[3] = 0;
  VAL[1] = sqrt(SMLNUM);
  VAL[2] = ONE;
  VAL[3] = sqrt(sqrt(BIGNUM));

  // Read input data until N=0.  Assume input eigenvalues are sorted
  // lexicographically (increasing by real part, then decreasing by
  // imaginary part)

  while (true) {
    (N, NDIM, ISRT) = await NIN.readInt3();
    if (N == 0) return;
    await NIN.readArray(ISELEC, NDIM);
    await NIN.readMatrix(TMP, N, N);
    (SIN, SEPIN) = await NIN.readDouble2();

    TNRM = zlange('M.value', N, N, TMP, LDT, RWORK);
    for (ISCL = 1; ISCL <= 3; ISCL++) {
      // 200

      // Scale input matrix

      KNT.value = KNT.value + 1;
      zlacpy('F', N, N, TMP, LDT, T, LDT);
      VMUL = VAL[ISCL];
      for (I = 1; I <= N; I++) {
        // 30
        zdscal(N, VMUL, T(1, I).asArray(), 1);
      } // 30
      if (TNRM == ZERO) VMUL = ONE;
      zlacpy('F', N, N, T, LDT, TSAV, LDT);

      // Compute Schur form

      zgehrd(N, 1, N, T, LDT, WORK(1), WORK(N + 1), LWORK - N, INFO);
      if (INFO.value != 0) {
        LMAX[1] = KNT.value;
        NINFO[1] = NINFO[1] + 1;
        continue;
      }

      // Generate unitary matrix

      zlacpy('L', N, N, T, LDT, Q, LDT);
      zunghr(N, 1, N, Q, LDT, WORK(1), WORK(N + 1), LWORK - N, INFO);

      // Compute Schur form

      for (J = 1; J <= N - 2; J++) {
        // 50
        for (I = J + 2; I <= N; I++) {
          // 40
          T[I][J] = Complex.zero;
        } // 40
      } // 50
      zhseqr('S.value', 'V', N, 1, N, T, LDT, W, Q, LDT, WORK, LWORK, INFO);
      if (INFO.value != 0) {
        LMAX[2] = KNT.value;
        NINFO[2] = NINFO[2] + 1;
        continue;
      }

      // Sort, select eigenvalues

      for (I = 1; I <= N; I++) {
        // 60
        IPNT[I] = I;
        SELECT[I] = false;
      } // 60
      if (ISRT == 0) {
        for (I = 1; I <= N; I++) {
          // 70
          WSRT[I] = (W[I]).toDouble();
        } // 70
      } else {
        for (I = 1; I <= N; I++) {
          // 80
          WSRT[I] = W[I].imaginary;
        } // 80
      }
      for (I = 1; I <= N - 1; I++) {
        // 100
        KMIN = I;
        VMIN = WSRT[I];
        for (J = I + 1; J <= N; J++) {
          // 90
          if (WSRT[J] < VMIN) {
            KMIN = J;
            VMIN = WSRT[J];
          }
        } // 90
        WSRT[KMIN] = WSRT[I];
        WSRT[I] = VMIN;
        ITMP = IPNT[I];
        IPNT[I] = IPNT[KMIN];
        IPNT[KMIN] = ITMP;
      } // 100
      for (I = 1; I <= NDIM; I++) {
        // 110
        SELECT[IPNT[ISELEC[I]]] = true;
      } // 110

      // Compute condition numbers

      zlacpy('F', N, N, Q, LDT, QSAV, LDT);
      zlacpy('F', N, N, T, LDT, TSAV1, LDT);
      ztrsen('B', 'V', SELECT, N, T, LDT, Q, LDT, WTMP, M, S, SEP, WORK, LWORK,
          INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      SEPTMP.value = SEP.value / VMUL;
      STMP.value = S.value;

      // Compute residuals

      zhst01(N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT);
      VMAX = max(RESULT[1], RESULT[2]);
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }

      // Compare condition number for eigenvalue cluster
      // taking its condition number into account

      V = max(TWO * N.toDouble() * EPS * TNRM, SMLNUM);
      if (TNRM == ZERO) V = ONE;
      if (V > SEPTMP.value) {
        TOL = ONE;
      } else {
        TOL = V / SEPTMP.value;
      }
      if (V > SEPIN) {
        TOLIN = ONE;
      } else {
        TOLIN = V / SEPIN;
      }
      TOL = max(TOL, SMLNUM / EPS);
      TOLIN = max(TOLIN, SMLNUM / EPS);
      if (EPS * (SIN - TOLIN) > STMP.value + TOL) {
        VMAX = ONE / EPS;
      } else if (SIN - TOLIN > STMP.value + TOL) {
        VMAX = (SIN - TOLIN) / (STMP.value + TOL);
      } else if (SIN + TOLIN < EPS * (STMP.value - TOL)) {
        VMAX = ONE / EPS;
      } else if (SIN + TOLIN < STMP.value - TOL) {
        VMAX = (STMP.value - TOL) / (SIN + TOLIN);
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[2]) {
        RMAX[2] = VMAX;
        if (NINFO[2] == 0) LMAX[2] = KNT.value;
      }

      // Compare condition numbers for invariant subspace
      // taking its condition number into account

      if (V > SEPTMP.value * STMP.value) {
        TOL = SEPTMP.value;
      } else {
        TOL = V / STMP.value;
      }
      if (V > SEPIN * SIN) {
        TOLIN = SEPIN;
      } else {
        TOLIN = V / SIN;
      }
      TOL = max(TOL, SMLNUM / EPS);
      TOLIN = max(TOLIN, SMLNUM / EPS);
      if (EPS * (SEPIN - TOLIN) > SEPTMP.value + TOL) {
        VMAX = ONE / EPS;
      } else if (SEPIN - TOLIN > SEPTMP.value + TOL) {
        VMAX = (SEPIN - TOLIN) / (SEPTMP.value + TOL);
      } else if (SEPIN + TOLIN < EPS * (SEPTMP.value - TOL)) {
        VMAX = ONE / EPS;
      } else if (SEPIN + TOLIN < SEPTMP.value - TOL) {
        VMAX = (SEPTMP.value - TOL) / (SEPIN + TOLIN);
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[2]) {
        RMAX[2] = VMAX;
        if (NINFO[2] == 0) LMAX[2] = KNT.value;
      }

      // Compare condition number for eigenvalue cluster
      // without taking its condition number into account

      if (SIN <= 2 * N * EPS && STMP.value <= 2 * N * EPS) {
        VMAX = ONE;
      } else if (EPS * SIN > STMP.value) {
        VMAX = ONE / EPS;
      } else if (SIN > STMP.value) {
        VMAX = SIN / STMP.value;
      } else if (SIN < EPS * STMP.value) {
        VMAX = ONE / EPS;
      } else if (SIN < STMP.value) {
        VMAX = STMP.value / SIN;
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[3]) {
        RMAX[3] = VMAX;
        if (NINFO[3] == 0) LMAX[3] = KNT.value;
      }

      // Compare condition numbers for invariant subspace
      // without taking its condition number into account

      if (SEPIN <= V && SEPTMP.value <= V) {
        VMAX = ONE;
      } else if (EPS * SEPIN > SEPTMP.value) {
        VMAX = ONE / EPS;
      } else if (SEPIN > SEPTMP.value) {
        VMAX = SEPIN / SEPTMP.value;
      } else if (SEPIN < EPS * SEPTMP.value) {
        VMAX = ONE / EPS;
      } else if (SEPIN < SEPTMP.value) {
        VMAX = SEPTMP.value / SEPIN;
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[3]) {
        RMAX[3] = VMAX;
        if (NINFO[3] == 0) LMAX[3] = KNT.value;
      }

      // Compute eigenvalue condition number only and compare
      // Update Q

      VMAX = ZERO;
      zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      zlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP.value = -ONE;
      STMP.value = -ONE;
      ztrsen('E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP,
          WORK, LWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (S.value != STMP.value) VMAX = ONE / EPS;
      if (-ONE != SEPTMP.value) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        // 130
        for (J = 1; J <= N; J++) {
          // 120
          if (TTMP(I, J) != T(I, J)) VMAX = ONE / EPS;
          if (QTMP(I, J) != Q(I, J)) VMAX = ONE / EPS;
        } // 120
      } // 130

      // Compute invariant subspace condition number only and compare
      // Update Q

      zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      zlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP.value = -ONE;
      STMP.value = -ONE;
      ztrsen('V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP,
          WORK, LWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (-ONE != STMP.value) VMAX = ONE / EPS;
      if (SEP.value != SEPTMP.value) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        // 150
        for (J = 1; J <= N; J++) {
          // 140
          if (TTMP(I, J) != T(I, J)) VMAX = ONE / EPS;
          if (QTMP(I, J) != Q(I, J)) VMAX = ONE / EPS;
        } // 140
      } // 150

      // Compute eigenvalue condition number only and compare
      // Do not update Q

      zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      zlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP.value = -ONE;
      STMP.value = -ONE;
      ztrsen('E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP,
          WORK, LWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (S.value != STMP.value) VMAX = ONE / EPS;
      if (-ONE != SEPTMP.value) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        // 170
        for (J = 1; J <= N; J++) {
          // 160
          if (TTMP(I, J) != T(I, J)) VMAX = ONE / EPS;
          if (QTMP(I, J) != QSAV(I, J)) VMAX = ONE / EPS;
        } // 160
      } // 170

      // Compute invariant subspace condition number only and compare
      // Do not update Q

      zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      zlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP.value = -ONE;
      STMP.value = -ONE;
      ztrsen('V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP,
          WORK, LWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (-ONE != STMP.value) VMAX = ONE / EPS;
      if (SEP.value != SEPTMP.value) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        // 190
        for (J = 1; J <= N; J++) {
          // 180
          if (TTMP(I, J) != T(I, J)) VMAX = ONE / EPS;
          if (QTMP(I, J) != QSAV(I, J)) VMAX = ONE / EPS;
        } // 180
      } // 190
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }
    } // 200
  }
}
