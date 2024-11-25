// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgehrd.dart';
import 'package:dart_lapack/src/dhseqr.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dorghr.dart';
import 'package:dart_lapack/src/dtrsen.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'dhst01.dart';

Future<void> dget38(
  final Array<double> RMAX_,
  final Array<int> LMAX_,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final Nin NIN,
  final TestDriver test,
  final String group,
  final double THRESH,
) async {
  final RMAX = RMAX_.having();
  final LMAX = LMAX_.having();
  final NINFO = NINFO_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8; // 2**(-24) = precision to which input data computed
  const LDT = 20, LWORK = 2 * LDT * (10 + LDT);
  const LIWORK = LDT * LDT;
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
  final WI = Array<double>(LDT),
      WITMP = Array<double>(LDT),
      WORK = Array<double>(LWORK),
      WR = Array<double>(LDT),
      RESULT = Array<double>(2),
      WRTMP = Array<double>(LDT);

  final SMLNUM = dlamch('S') / dlamch('P');
  final BIGNUM = ONE / SMLNUM;
  final EPS = max(dlamch('P'), EPSIN);

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

  final VAL = [sqrt(SMLNUM), ONE, sqrt(sqrt(BIGNUM))];

  // Read input data until N=0.  Assume input eigenvalues are sorted
  // lexicographically (increasing by real part, then decreasing by
  // imaginary part)

  while (true) {
    final (N, NDIM) = await NIN.readInt2();
    if (N == 0) return;
    await NIN.readArray(ISELEC, NDIM);
    await NIN.readMatrix(TMP, N, N);
    final (SIN, SEPIN) = await NIN.readDouble2();
    final TNRM = dlange('M', N, N, TMP, LDT, WORK);

    final ctx = (TMP: TMP.copy(), ISELEC: ISELEC.copy());
    test.group(group, () {
      test('DGET38', () {
        final (:TMP, :ISELEC) = ctx;
        final INFO = Box(0), M = Box(0);

        for (var ISCL = 1; ISCL <= 3; ISCL++) {
          // Scale input matrix

          KNT.value++;
          dlacpy('F', N, N, TMP, LDT, T, LDT);
          var VMUL = VAL[ISCL - 1];
          for (var I = 1; I <= N; I++) {
            dscal(N, VMUL, T(1, I).asArray(), 1);
          }
          if (TNRM == ZERO) VMUL = ONE;
          dlacpy('F', N, N, T, LDT, TSAV, LDT);

          // Compute Schur form

          dgehrd(N, 1, N, T, LDT, WORK, WORK(N + 1), LWORK - N, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[1] = KNT.value;
            NINFO[1]++;
            continue;
          }

          // Generate orthogonal matrix

          dlacpy('L', N, N, T, LDT, Q, LDT);
          dorghr(N, 1, N, Q, LDT, WORK(1), WORK(N + 1), LWORK - N, INFO);

          // Compute Schur form

          dhseqr('S', 'V', N, 1, N, T, LDT, WR, WI, Q, LDT, WORK, LWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[2] = KNT.value;
            NINFO[2]++;
            continue;
          }

          // Sort, select eigenvalues

          for (var I = 1; I <= N; I++) {
            IPNT[I] = I;
            SELECT[I] = false;
          }
          dcopy(N, WR, 1, WRTMP, 1);
          dcopy(N, WI, 1, WITMP, 1);
          for (var I = 1; I <= N - 1; I++) {
            var KMIN = I;
            var VRMIN = WRTMP[I];
            var VIMIN = WITMP[I];
            for (var J = I + 1; J <= N; J++) {
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
            final ITMP = IPNT[I];
            IPNT[I] = IPNT[KMIN];
            IPNT[KMIN] = ITMP;
          }
          for (var I = 1; I <= NDIM; I++) {
            SELECT[IPNT[ISELEC[I]]] = true;
          }

          // Compute condition numbers

          dlacpy('F', N, N, Q, LDT, QSAV, LDT);
          dlacpy('F', N, N, T, LDT, TSAV1, LDT);
          final S = Box(0.0),
              SEP = Box(0.0),
              SEPTMP = Box(0.0),
              STMP = Box(0.0);
          dtrsen('B', 'V', SELECT, N, T, LDT, Q, LDT, WRTMP, WITMP, M, S, SEP,
              WORK, LWORK, IWORK, LIWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          SEPTMP.value = SEP.value / VMUL;
          STMP.value = S.value;

          // Compute residuals
          {
            dhst01(N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RESULT);
            final VMAX = max(RESULT[1], RESULT[2]);
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[1]) {
              RMAX[1] = VMAX;
              if (NINFO[1] == 0) LMAX[1] = KNT.value;
            }
          }

          final V = TNRM == ZERO ? ONE : max(TWO * N * EPS * TNRM, SMLNUM);

          // Compare condition number for eigenvalue cluster
          // taking its condition number into account
          {
            double TOL;
            if (V > SEPTMP.value) {
              TOL = ONE;
            } else {
              TOL = V / SEPTMP.value;
            }
            double TOLIN;
            if (V > SEPIN) {
              TOLIN = ONE;
            } else {
              TOLIN = V / SEPIN;
            }
            TOL = max(TOL, SMLNUM / EPS);
            TOLIN = max(TOLIN, SMLNUM / EPS);
            final double VMAX;
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
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[2]) {
              RMAX[2] = VMAX;
              if (NINFO[2] == 0) LMAX[2] = KNT.value;
            }
          }

          // Compare condition numbers for invariant subspace
          // taking its condition number into account
          {
            double TOL;
            if (V > SEPTMP.value * STMP.value) {
              TOL = SEPTMP.value;
            } else {
              TOL = V / STMP.value;
            }
            double TOLIN;
            if (V > SEPIN * SIN) {
              TOLIN = SEPIN;
            } else {
              TOLIN = V / SIN;
            }
            TOL = max(TOL, SMLNUM / EPS);
            TOLIN = max(TOLIN, SMLNUM / EPS);
            final double VMAX;
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
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[2]) {
              RMAX[2] = VMAX;
              if (NINFO[2] == 0) LMAX[2] = KNT.value;
            }
          }

          // Compare condition number for eigenvalue cluster
          // without taking its condition number into account
          {
            final double VMAX;
            if (SIN <= (2 * N) * EPS && STMP.value <= (2 * N) * EPS) {
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
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[3]) {
              RMAX[3] = VMAX;
              if (NINFO[3] == 0) LMAX[3] = KNT.value;
            }
          }

          // Compare condition numbers for invariant subspace
          // without taking its condition number into account
          {
            final double VMAX;
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
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[3]) {
              RMAX[3] = VMAX;
              if (NINFO[3] == 0) LMAX[3] = KNT.value;
            }
          }

          // Compute eigenvalue condition number only and compare
          // Update Q

          var VMAX = ZERO;
          dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
          dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
          SEPTMP.value = -ONE;
          STMP.value = -ONE;
          dtrsen('E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M,
              STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          if (S.value != STMP.value) VMAX = ONE / EPS;
          if (-ONE != SEPTMP.value) VMAX = ONE / EPS;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
              if (QTMP[I][J] != Q[I][J]) VMAX = ONE / EPS;
            }
          }

          // Compute invariant subspace condition number only and compare
          // Update Q

          dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
          dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
          SEPTMP.value = -ONE;
          STMP.value = -ONE;
          dtrsen('V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M,
              STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          if (-ONE != STMP.value) VMAX = ONE / EPS;
          if (SEP.value != SEPTMP.value) VMAX = ONE / EPS;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
              if (QTMP[I][J] != Q[I][J]) VMAX = ONE / EPS;
            }
          }

          // Compute eigenvalue condition number only and compare
          // Do not update Q

          dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
          dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
          SEPTMP.value = -ONE;
          STMP.value = -ONE;
          dtrsen('E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M,
              STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          if (S.value != STMP.value) VMAX = ONE / EPS;
          if (-ONE != SEPTMP.value) VMAX = ONE / EPS;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
              if (QTMP[I][J] != QSAV[I][J]) VMAX = ONE / EPS;
            }
          }

          // Compute invariant subspace condition number only and compare
          // Do not update Q

          dlacpy('F', N, N, TSAV1, LDT, TTMP, LDT);
          dlacpy('F', N, N, QSAV, LDT, QTMP, LDT);
          SEPTMP.value = -ONE;
          STMP.value = -ONE;
          dtrsen('V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M,
              STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          if (-ONE != STMP.value) VMAX = ONE / EPS;
          if (SEP.value != SEPTMP.value) VMAX = ONE / EPS;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              if (TTMP[I][J] != T[I][J]) VMAX = ONE / EPS;
              if (QTMP[I][J] != QSAV[I][J]) VMAX = ONE / EPS;
            }
          }
          test.expect(VMAX, lessThan(THRESH));
          if (VMAX > RMAX[1]) {
            RMAX[1] = VMAX;
            if (NINFO[1] == 0) LMAX[1] = KNT.value;
          }
        }
      });
    });
  }
}
