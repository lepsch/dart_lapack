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
import 'package:dart_lapack/src/dtrevc.dart';
import 'package:dart_lapack/src/dtrsna.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';

Future<void> dget37(
  final Array<double> RMAX_,
  final Array<int> LMAX_,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final Nin NIN,
  final TestDriver test,
  final String group,
  final double THRESH,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RMAX = RMAX_.having();
  final LMAX = LMAX_.having();
  final NINFO = NINFO_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8; // 2**(-24) = precision to which input data computed
  const LDT = 20, LWORK = 2 * LDT * (10 + LDT);
  final SELECT = Array<bool>(LDT);
  final IWORK = Array<int>(2 * LDT), LCMP = Array<int>(3);
  final LE = Matrix<double>(LDT, LDT),
      RE = Matrix<double>(LDT, LDT),
      T = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT);

  var DUM = Array<double>(1),
      S = Array<double>(LDT),
      SEP = Array<double>(LDT),
      SEPIN = Array<double>(LDT),
      SEPTMP = Array<double>(LDT),
      SIN = Array<double>(LDT),
      STMP = Array<double>(LDT),
      WI = Array<double>(LDT),
      WITMP = Array<double>(LDT),
      WORK = Array<double>(LWORK),
      WR = Array<double>(LDT),
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

  final VAL = [sqrt(SMLNUM), ONE, sqrt(BIGNUM)];

  // Read input data until N=0.  Assume input eigenvalues are sorted
  // lexicographically (increasing by real part, then decreasing by
  // imaginary part)

  while (true) {
    final N = await NIN.readInt();
    if (N == 0) return;
    await NIN.readMatrix(TMP, N, N);
    for (var I = 1; I <= N; I++) {
      final (_, _, f3, f4) = await NIN.readDouble4();
      SIN[I] = f3;
      SEPIN[I] = f4;
    }
    final TNRM = dlange('M', N, N, TMP, LDT, WORK);

    final ctx = (TMP: TMP.copy(), SIN: SIN.copy(), SEPIN: SEPIN.copy());

    test.group(group, () {
      test('DGET37', () {
        final (:TMP, :SIN, :SEPIN) = ctx;
        final INFO = Box(0), M = Box(0);

        // Begin test
        for (var ISCL = 1; ISCL <= 3; ISCL++) {
          // Scale input matrix
          KNT.value++;
          dlacpy('F', N, N, TMP, LDT, T, LDT);
          var VMUL = VAL[ISCL - 1];
          for (var I = 1; I <= N; I++) {
            dscal(N, VMUL, T(1, I).asArray(), 1);
          }
          if (TNRM == ZERO) VMUL = ONE;

          // Compute eigenvalues and eigenvectors
          dgehrd(N, 1, N, T, LDT, WORK(1), WORK(N + 1), LWORK - N, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[1] = KNT.value;
            NINFO[1]++;
            continue;
          }
          for (var J = 1; J <= N - 2; J++) {
            for (var I = J + 2; I <= N; I++) {
              T[I][J] = ZERO;
            }
          }

          // Compute Schur form

          dhseqr('S', 'N', N, 1, N, T, LDT, WR, WI, DUM.asMatrix(1), 1, WORK,
              LWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[2] = KNT.value;
            NINFO[2]++;
            continue;
          }

          // Compute eigenvectors

          dtrevc('Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, N, M, WORK,
              INFO);

          // Compute condition numbers

          dtrsna('Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, S, SEP, N,
              M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }

          // Sort eigenvalues and condition numbers lexicographically
          // to compare with inputs

          dcopy(N, WR, 1, WRTMP, 1);
          dcopy(N, WI, 1, WITMP, 1);
          dcopy(N, S, 1, STMP, 1);
          dcopy(N, SEP, 1, SEPTMP, 1);
          dscal(N, ONE / VMUL, SEPTMP, 1);
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
            VRMIN = STMP[KMIN];
            STMP[KMIN] = STMP[I];
            STMP[I] = VRMIN;
            VRMIN = SEPTMP[KMIN];
            SEPTMP[KMIN] = SEPTMP[I];
            SEPTMP[I] = VRMIN;
          }

          // Compare condition numbers for eigenvalues
          // taking their condition numbers into account

          final V = TNRM == ZERO ? ONE : max(TWO * N * EPS * TNRM, SMLNUM);
          for (var I = 1; I <= N; I++) {
            double TOL;
            if (V > SEPTMP[I]) {
              TOL = ONE;
            } else {
              TOL = V / SEPTMP[I];
            }
            double TOLIN;
            if (V > SEPIN[I]) {
              TOLIN = ONE;
            } else {
              TOLIN = V / SEPIN[I];
            }
            TOL = max(TOL, SMLNUM / EPS);
            TOLIN = max(TOLIN, SMLNUM / EPS);
            final double VMAX;
            if (EPS * (SIN[I] - TOLIN) > STMP[I] + TOL) {
              VMAX = ONE / EPS;
            } else if (SIN[I] - TOLIN > STMP[I] + TOL) {
              VMAX = (SIN[I] - TOLIN) / (STMP[I] + TOL);
            } else if (SIN[I] + TOLIN < EPS * (STMP[I] - TOL)) {
              VMAX = ONE / EPS;
            } else if (SIN[I] + TOLIN < STMP[I] - TOL) {
              VMAX = (STMP[I] - TOL) / (SIN[I] + TOLIN);
            } else {
              VMAX = ONE;
            }
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[2]) {
              RMAX[2] = VMAX;
              if (NINFO[2] == 0) LMAX[2] = KNT.value;
            }
          }

          // Compare condition numbers for eigenvectors
          // taking their condition numbers into account

          for (var I = 1; I <= N; I++) {
            double TOL;
            if (V > SEPTMP[I] * STMP[I]) {
              TOL = SEPTMP[I];
            } else {
              TOL = V / STMP[I];
            }
            double TOLIN;
            if (V > SEPIN[I] * SIN[I]) {
              TOLIN = SEPIN[I];
            } else {
              TOLIN = V / SIN[I];
            }
            TOL = max(TOL, SMLNUM / EPS);
            TOLIN = max(TOLIN, SMLNUM / EPS);
            final double VMAX;
            if (EPS * (SEPIN[I] - TOLIN) > SEPTMP[I] + TOL) {
              VMAX = ONE / EPS;
            } else if (SEPIN[I] - TOLIN > SEPTMP[I] + TOL) {
              VMAX = (SEPIN[I] - TOLIN) / (SEPTMP[I] + TOL);
            } else if (SEPIN[I] + TOLIN < EPS * (SEPTMP[I] - TOL)) {
              VMAX = ONE / EPS;
            } else if (SEPIN[I] + TOLIN < SEPTMP[I] - TOL) {
              VMAX = (SEPTMP[I] - TOL) / (SEPIN[I] + TOLIN);
            } else {
              VMAX = ONE;
            }
            test.expect(VMAX, lessThan(THRESH));
            if (VMAX > RMAX[2]) {
              RMAX[2] = VMAX;
              if (NINFO[2] == 0) LMAX[2] = KNT.value;
            }
          }

          // Compare condition numbers for eigenvalues
          // without taking their condition numbers into account

          for (var I = 1; I <= N; I++) {
            final double VMAX;
            if (SIN[I] <= (2 * N) * EPS && STMP[I] <= (2 * N) * EPS) {
              VMAX = ONE;
            } else if (EPS * SIN[I] > STMP[I]) {
              VMAX = ONE / EPS;
            } else if (SIN[I] > STMP[I]) {
              VMAX = SIN[I] / STMP[I];
            } else if (SIN[I] < EPS * STMP[I]) {
              VMAX = ONE / EPS;
            } else if (SIN[I] < STMP[I]) {
              VMAX = STMP[I] / SIN[I];
            } else {
              VMAX = ONE;
            }
            if (VMAX > RMAX[3]) {
              RMAX[3] = VMAX;
              if (NINFO[3] == 0) LMAX[3] = KNT.value;
            }
          }

          // Compare condition numbers for eigenvectors
          // without taking their condition numbers into account

          for (var I = 1; I <= N; I++) {
            final double VMAX;
            if (SEPIN[I] <= V && SEPTMP[I] <= V) {
              VMAX = ONE;
            } else if (EPS * SEPIN[I] > SEPTMP[I]) {
              VMAX = ONE / EPS;
            } else if (SEPIN[I] > SEPTMP[I]) {
              VMAX = SEPIN[I] / SEPTMP[I];
            } else if (SEPIN[I] < EPS * SEPTMP[I]) {
              VMAX = ONE / EPS;
            } else if (SEPIN[I] < SEPTMP[I]) {
              VMAX = SEPTMP[I] / SEPIN[I];
            } else {
              VMAX = ONE;
            }
            if (VMAX > RMAX[3]) {
              RMAX[3] = VMAX;
              if (NINFO[3] == 0) LMAX[3] = KNT.value;
            }
          }

          // Compute eigenvalue condition numbers only and compare

          var VMAX = ZERO;
          DUM[1] = -ONE;
          dcopy(N, DUM, 0, STMP, 1);
          dcopy(N, DUM, 0, SEPTMP, 1);
          dtrsna('Eigcond', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= N; I++) {
            if (STMP[I] != S[I]) VMAX = ONE / EPS;
            if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
          }

          // Compute eigenvector condition numbers only and compare

          dcopy(N, DUM, 0, STMP, 1);
          dcopy(N, DUM, 0, SEPTMP, 1);
          dtrsna('Veccond', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= N; I++) {
            if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
            if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
          }

          // Compute all condition numbers using SELECT and compare

          for (var I = 1; I <= N; I++) {
            SELECT[I] = true;
          }
          dcopy(N, DUM, 0, STMP, 1);
          dcopy(N, DUM, 0, SEPTMP, 1);
          dtrsna('Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= N; I++) {
            if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
            if (STMP[I] != S[I]) VMAX = ONE / EPS;
          }

          // Compute eigenvalue condition numbers using SELECT and compare

          dcopy(N, DUM, 0, STMP, 1);
          dcopy(N, DUM, 0, SEPTMP, 1);
          dtrsna('Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= N; I++) {
            if (STMP[I] != S[I]) VMAX = ONE / EPS;
            if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
          }

          // Compute eigenvector condition numbers using SELECT and compare

          dcopy(N, DUM, 0, STMP, 1);
          dcopy(N, DUM, 0, SEPTMP, 1);
          dtrsna('Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= N; I++) {
            if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
            if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
          }
          test.expect(VMAX, lessThan(THRESH));
          if (VMAX > RMAX[1]) {
            RMAX[1] = VMAX;
            if (NINFO[1] == 0) LMAX[1] = KNT.value;
          }

          // Select first real and first complex eigenvalue

          final int ICMP;
          if (WI[1] == ZERO) {
            LCMP[1] = 1;
            var IFND = 0;
            for (var I = 2; I <= N; I++) {
              if (IFND == 1 || WI[I] == ZERO) {
                SELECT[I] = false;
              } else {
                IFND = 1;
                LCMP[2] = I;
                LCMP[3] = I + 1;
                dcopy(N, RE(1, I).asArray(), 1, RE(1, 2).asArray(), 1);
                dcopy(N, RE(1, I + 1).asArray(), 1, RE(1, 3).asArray(), 1);
                dcopy(N, LE(1, I).asArray(), 1, LE(1, 2).asArray(), 1);
                dcopy(N, LE(1, I + 1).asArray(), 1, LE(1, 3).asArray(), 1);
              }
            }
            if (IFND == 0) {
              ICMP = 1;
            } else {
              ICMP = 3;
            }
          } else {
            LCMP[1] = 1;
            LCMP[2] = 2;
            var IFND = 0;
            for (var I = 3; I <= N; I++) {
              if (IFND == 1 || WI[I] != ZERO) {
                SELECT[I] = false;
              } else {
                LCMP[3] = I;
                IFND = 1;
                dcopy(N, RE(1, I).asArray(), 1, RE(1, 3).asArray(), 1);
                dcopy(N, LE(1, I).asArray(), 1, LE(1, 3).asArray(), 1);
              }
            }
            if (IFND == 0) {
              ICMP = 2;
            } else {
              ICMP = 3;
            }
          }

          // Compute all selected condition numbers

          dcopy(ICMP, DUM, 0, STMP, 1);
          dcopy(ICMP, DUM, 0, SEPTMP, 1);
          dtrsna('Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= ICMP; I++) {
            final J = LCMP[I];
            if (SEPTMP[I] != SEP[J]) VMAX = ONE / EPS;
            if (STMP[I] != S[J]) VMAX = ONE / EPS;
          }

          // Compute selected eigenvalue condition numbers

          dcopy(ICMP, DUM, 0, STMP, 1);
          dcopy(ICMP, DUM, 0, SEPTMP, 1);
          dtrsna('Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= ICMP; I++) {
            final J = LCMP[I];
            if (STMP[I] != S[J]) VMAX = ONE / EPS;
            if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
          }

          // Compute selected eigenvector condition numbers

          dcopy(ICMP, DUM, 0, STMP, 1);
          dcopy(ICMP, DUM, 0, SEPTMP, 1);
          dtrsna('Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP,
              SEPTMP, N, M, WORK.asMatrix(N), N, IWORK, INFO);
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            LMAX[3] = KNT.value;
            NINFO[3]++;
            continue;
          }
          for (var I = 1; I <= ICMP; I++) {
            final J = LCMP[I];
            if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
            if (SEPTMP[I] != SEP[J]) VMAX = ONE / EPS;
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
