import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgehrd.dart';
import 'package:lapack/src/dhseqr.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dorghr.dart';
import 'package:lapack/src/dtrsen.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'dhst01.dart';

Future<void> dget38(
  final Array<double> RMAX,
  final Array<int> LMAX,
  final Array<int> NINFO,
  final Box<int> KNT,
  final Nin NIN,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8;
  const LDT = 20, LWORK = 2 * LDT * (10 + LDT);
  const LIWORK = LDT * LDT;
  int I, ISCL, ITMP, J, KMIN, M = 0, N = 0, NDIM = 0;
  double BIGNUM,
      EPS,
      S = 0,
      SEP = 0,
      SEPIN = 0,
      SEPTMP,
      SIN = 0,
      SMLNUM,
      STMP,
      TNRM,
      TOL,
      TOLIN,
      V,
      VIMIN,
      VMAX,
      VMUL,
      VRMIN;
  final SELECT = Array<bool>(LDT);
  final IPNT = Array<int>(LDT),
      ISELEC = Array<int>(LDT),
      IWORK = Array<int>(LIWORK);
  final Q = Matrix<double>(LDT, LDT),
      QSAV = Matrix<double>(LDT, LDT),
      QTMP = Matrix<double>(LDT, LDT),
      T = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT),
      TSAV = Matrix<double>(LDT, LDT),
      TSAV1 = Matrix<double>(LDT, LDT),
      TTMP = Matrix<double>(LDT, LDT);

  final VAL = Array<double>(3),
      WI = Array<double>(LDT),
      WITMP = Array<double>(LDT),
      WORK = Array<double>(LWORK),
      WR = Array<double>(LDT),
      RESULT = Array<double>(2),
      WRTMP = Array<double>(LDT);

  final INFO = Box(0);

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
    (N, NDIM) = await NIN.readInt2();
    if (N == 0) return;
    await NIN.readArray(ISELEC, NDIM);
    await NIN.readMatrix(TMP, N, N);
    for (I = 1; I <= N; I++) {}
    (SIN, SEPIN) = await NIN.readDouble2();

    TNRM = dlange('M', N, N, TMP, LDT, WORK);
    for (ISCL = 1; ISCL <= 3; ISCL++) {
      // Scale input matrix

      KNT.value = KNT.value + 1;
      dlacpy('F', N, N, TMP, LDT, T, LDT);
      VMUL = VAL[ISCL];
      for (I = 1; I <= N; I++) {
        dscal(N, VMUL, T(1, I).asArray(), 1);
      }
      if (TNRM == ZERO) VMUL = ONE;
      dlacpy('F', N, N, T, LDT, TSAV, LDT);

      // Compute Schur form

      dgehrd(N, 1, N, T, LDT, WORK, WORK(N + 1), LWORK - N, INFO);
      if (INFO.value != 0) {
        LMAX[1] = KNT.value;
        NINFO[1] = NINFO[1] + 1;
        continue;
      }

      // Generate orthogonal matrix

      dlacpy('L', N, N, T, LDT, Q, LDT);
      dorghr(N, 1, N, Q, LDT, WORK[1], WORK[N + 1], LWORK - N, INFO.value);

      // Compute Schur form

      dhseqr(
        'S',
        'V',
        N,
        1,
        N,
        T,
        LDT,
        WR,
        WI,
        Q,
        LDT,
        WORK,
        LWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[2] = KNT.value;
        NINFO[2] = NINFO[2] + 1;
        continue;
      }

      // Sort, select eigenvalues

      for (I = 1; I <= N; I++) {
        IPNT[I] = I;
        SELECT[I] = false;
      }
      dcopy(N, WR, 1, WRTMP, 1);
      dcopy(N, WI, 1, WITMP, 1);
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
      for (I = 1; I <= NDIM; I++) {
        SELECT[IPNT[ISELEC[I]]] = true;
      }

      // Compute condition numbers

      dlacpy('F', N, N, Q, LDT, QSAV, LDT);
      dlacpy('F', N, N, T, LDT, TSAV1, LDT);
      dtrsen(
        'B',
        'V',
        SELECT,
        N,
        T,
        LDT,
        Q,
        LDT,
        WRTMP,
        WITMP,
        M,
        S,
        SEP,
        WORK,
        LWORK,
        IWORK,
        LIWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      SEPTMP = SEP / VMUL;
      STMP = S;

      // Compute residuals

      dhst01(N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RESULT);
      VMAX = max(RESULT[1], RESULT[2]);
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }

      // Compare condition number for eigenvalue cluster
      // taking its condition number into account

      V = max(TWO * N.toDouble() * EPS * TNRM, SMLNUM);
      if (TNRM == ZERO) V = ONE;
      if (V > SEPTMP) {
        TOL = ONE;
      } else {
        TOL = V / SEPTMP;
      }
      if (V > SEPIN) {
        TOLIN = ONE;
      } else {
        TOLIN = V / SEPIN;
      }
      TOL = max(TOL, SMLNUM / EPS);
      TOLIN = max(TOLIN, SMLNUM / EPS);
      if (EPS * (SIN - TOLIN) > STMP + TOL) {
        VMAX = ONE / EPS;
      } else if (SIN - TOLIN > STMP + TOL) {
        VMAX = (SIN - TOLIN) / (STMP + TOL);
      } else if (SIN + TOLIN < EPS * (STMP - TOL)) {
        VMAX = ONE / EPS;
      } else if (SIN + TOLIN < STMP - TOL) {
        VMAX = (STMP - TOL) / (SIN + TOLIN);
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[2]) {
        RMAX[2] = VMAX;
        if (NINFO[2] == 0) LMAX[2] = KNT.value;
      }

      // Compare condition numbers for invariant subspace
      // taking its condition number into account

      if (V > SEPTMP * STMP) {
        TOL = SEPTMP;
      } else {
        TOL = V / STMP;
      }
      if (V > SEPIN * SIN) {
        TOLIN = SEPIN;
      } else {
        TOLIN = V / SIN;
      }
      TOL = max(TOL, SMLNUM / EPS);
      TOLIN = max(TOLIN, SMLNUM / EPS);
      if (EPS * (SEPIN - TOLIN) > SEPTMP + TOL) {
        VMAX = ONE / EPS;
      } else if (SEPIN - TOLIN > SEPTMP + TOL) {
        VMAX = (SEPIN - TOLIN) / (SEPTMP + TOL);
      } else if (SEPIN + TOLIN < EPS * (SEPTMP - TOL)) {
        VMAX = ONE / EPS;
      } else if (SEPIN + TOLIN < SEPTMP - TOL) {
        VMAX = (SEPTMP - TOL) / (SEPIN + TOLIN);
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[2]) {
        RMAX[2] = VMAX;
        if (NINFO[2] == 0) LMAX[2] = KNT.value;
      }

      // Compare condition number for eigenvalue cluster
      // without taking its condition number into account

      if (SIN <= (2 * N).toDouble() * EPS && STMP <= (2 * N).toDouble() * EPS) {
        VMAX = ONE;
      } else if (EPS * SIN > STMP) {
        VMAX = ONE / EPS;
      } else if (SIN > STMP) {
        VMAX = SIN / STMP;
      } else if (SIN < EPS * STMP) {
        VMAX = ONE / EPS;
      } else if (SIN < STMP) {
        VMAX = STMP / SIN;
      } else {
        VMAX = ONE;
      }
      if (VMAX > RMAX[3]) {
        RMAX[3] = VMAX;
        if (NINFO[3] == 0) LMAX[3] = KNT.value;
      }

      // Compare condition numbers for invariant subspace
      // without taking its condition number into account

      if (SEPIN <= V && SEPTMP <= V) {
        VMAX = ONE;
      } else if (EPS * SEPIN > SEPTMP) {
        VMAX = ONE / EPS;
      } else if (SEPIN > SEPTMP) {
        VMAX = SEPIN / SEPTMP;
      } else if (SEPIN < EPS * SEPTMP) {
        VMAX = ONE / EPS;
      } else if (SEPIN < SEPTMP) {
        VMAX = SEPTMP / SEPIN;
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
      dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP = -ONE;
      STMP = -ONE;
      dtrsen(
        'E',
        'V',
        SELECT,
        N,
        TTMP,
        LDT,
        QTMP,
        LDT,
        WRTMP,
        WITMP,
        M,
        STMP,
        SEPTMP,
        WORK,
        LWORK,
        IWORK,
        LIWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (S != STMP) VMAX = ONE / EPS;
      if (-ONE != SEPTMP) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        for (J = 1; J <= N; J++) {
          if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
          if (QTMP[I][J] != Q[I][J]) VMAX = ONE / EPS;
        }
      }

      // Compute invariant subspace condition number only and compare
      // Update Q

      dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP = -ONE;
      STMP = -ONE;
      dtrsen(
        'V',
        'V',
        SELECT,
        N,
        TTMP,
        LDT,
        QTMP,
        LDT,
        WRTMP,
        WITMP,
        M,
        STMP,
        SEPTMP,
        WORK,
        LWORK,
        IWORK,
        LIWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (-ONE != STMP) VMAX = ONE / EPS;
      if (SEP != SEPTMP) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        for (J = 1; J <= N; J++) {
          if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
          if (QTMP[I][J] != Q[I][J]) VMAX = ONE / EPS;
        }
      }

      // Compute eigenvalue condition number only and compare
      // Do not update Q

      dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP = -ONE;
      STMP = -ONE;
      dtrsen(
        'E',
        'N',
        SELECT,
        N,
        TTMP,
        LDT,
        QTMP,
        LDT,
        WRTMP,
        WITMP,
        M,
        STMP,
        SEPTMP,
        WORK,
        LWORK,
        IWORK,
        LIWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (S != STMP) VMAX = ONE / EPS;
      if (-ONE != SEPTMP) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        for (J = 1; J <= N; J++) {
          if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
          if (QTMP[I][J] != QSAV[I][J]) VMAX = ONE / EPS;
        }
      }

      // Compute invariant subspace condition number only and compare
      // Do not update Q

      dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
      dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
      SEPTMP = -ONE;
      STMP = -ONE;
      dtrsen(
        'V',
        'N',
        SELECT,
        N,
        TTMP,
        LDT,
        QTMP,
        LDT,
        WRTMP,
        WITMP,
        M,
        STMP,
        SEPTMP,
        WORK,
        LWORK,
        IWORK,
        LIWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      if (-ONE != STMP) VMAX = ONE / EPS;
      if (SEP != SEPTMP) VMAX = ONE / EPS;
      for (I = 1; I <= N; I++) {
        for (J = 1; J <= N; J++) {
          if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
          if (QTMP[I][J] != QSAV[I][J]) VMAX = ONE / EPS;
        }
      }
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }
    }
  }
}
