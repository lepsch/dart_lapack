// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgehrd.dart';
import 'package:dart_lapack/src/zhseqr.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/ztrevc.dart';
import 'package:dart_lapack/src/ztrsna.dart';

Future<void> zget37(
  final Array<double> RMAX_,
  final Array<int> LMAX_,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final Nin NIN,
) async {
  final RMAX = RMAX_.having(length: 3);
  final LMAX = LMAX_.having(length: 3);
  final NINFO = NINFO_.having(length: 3);
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8;
  const LDT = 20, LWORK = 2 * LDT * (10 + LDT);
  int I, ICMP, ISCL, ISRT, J, KMIN, N;
  double BIGNUM, EPS, SMLNUM, TNRM, TOL, TOLIN, V, VCMIN, VMAX, VMIN, VMUL;
  final SELECT = Array<bool>(LDT);
  final LCMP = Array<int>(3);
  final DUM = Array<double>(1),
      RWORK = Array<double>(2 * LDT),
      S = Array<double>(LDT),
      SEP = Array<double>(LDT),
      SEPIN = Array<double>(LDT),
      SEPTMP = Array<double>(LDT),
      SIN = Array<double>(LDT),
      STMP = Array<double>(LDT),
      VAL = Array<double>(3),
      WIIN = Array<double>(LDT),
      WRIN = Array<double>(LDT),
      WSRT = Array<double>(LDT);
  final CDUM = Array<Complex>(1),
      W = Array<Complex>(LDT),
      WORK = Array<Complex>(LWORK),
      WTMP = Array<Complex>(LDT);
  final LE = Matrix<Complex>(LDT, LDT),
      RE = Matrix<Complex>(LDT, LDT),
      T = Matrix<Complex>(LDT, LDT),
      TMP = Matrix<Complex>(LDT, LDT);
  final INFO = Box(0), M = Box(0);

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
  VAL[3] = sqrt(BIGNUM);

  // Read input data until N=0.  Assume input eigenvalues are sorted
  // lexicographically (increasing by real part if ISRT = 0,
  // increasing by imaginary part if ISRT = 1)

  while (true) {
    (N, ISRT) = await NIN.readInt2();
    if (N == 0) return;
    await NIN.readMatrix(TMP, N, N);

    for (I = 1; I <= N; I++) {
      await NIN.readBoxes(WRIN(I), WIIN(I), SIN(I), SEPIN(I));
    }
    TNRM = zlange('M', N, N, TMP, LDT, RWORK);
    for (ISCL = 1; ISCL <= 3; ISCL++) {
      // Scale input matrix

      KNT.value++;
      zlacpy('F', N, N, TMP, LDT, T, LDT);
      VMUL = VAL[ISCL];
      for (I = 1; I <= N; I++) {
        zdscal(N, VMUL, T(1, I).asArray(), 1);
      }
      if (TNRM == ZERO) VMUL = ONE;

      // Compute eigenvalues and eigenvectors

      zgehrd(N, 1, N, T, LDT, WORK(1), WORK(N + 1), LWORK - N, INFO);
      if (INFO.value != 0) {
        LMAX[1] = KNT.value;
        NINFO[1]++;
        continue;
      }
      for (J = 1; J <= N - 2; J++) {
        for (I = J + 2; I <= N; I++) {
          T[I][J] = Complex.zero;
        }
      }

      // Compute Schur form

      zhseqr(
          'S', 'N', N, 1, N, T, LDT, W, CDUM.asMatrix(), 1, WORK, LWORK, INFO);
      if (INFO.value != 0) {
        LMAX[2] = KNT.value;
        NINFO[2]++;
        continue;
      }

      // Compute eigenvectors

      for (I = 1; I <= N; I++) {
        SELECT[I] = true;
      }
      ztrevc('B', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, N, M, WORK, RWORK,
          INFO);

      // Compute condition numbers

      ztrsna('B', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, S, SEP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }

      // Sort eigenvalues and condition numbers lexicographically
      // to compare with inputs

      zcopy(N, W, 1, WTMP, 1);
      if (ISRT == 0) {
        // Sort by increasing real part

        for (I = 1; I <= N; I++) {
          WSRT[I] = W[I].real;
        }
      } else {
        // Sort by increasing imaginary part

        for (I = 1; I <= N; I++) {
          WSRT[I] = W[I].imaginary;
        }
      }
      dcopy(N, S, 1, STMP, 1);
      dcopy(N, SEP, 1, SEPTMP, 1);
      dscal(N, ONE / VMUL, SEPTMP, 1);
      for (I = 1; I <= N - 1; I++) {
        KMIN = I;
        VMIN = WSRT[I];
        for (J = I + 1; J <= N; J++) {
          if (WSRT[J] < VMIN) {
            KMIN = J;
            VMIN = WSRT[J];
          }
        }
        WSRT[KMIN] = WSRT[I];
        WSRT[I] = VMIN;
        VCMIN = WTMP[I].real;
        WTMP[I] = W[KMIN];
        WTMP[KMIN] = VCMIN.toComplex();
        VMIN = STMP[KMIN];
        STMP[KMIN] = STMP[I];
        STMP[I] = VMIN;
        VMIN = SEPTMP[KMIN];
        SEPTMP[KMIN] = SEPTMP[I];
        SEPTMP[I] = VMIN;
      }

      // Compare condition numbers for eigenvalues
      // taking their condition numbers into account

      V = max(TWO * N * EPS * TNRM, SMLNUM);
      if (TNRM == ZERO) V = ONE;
      for (I = 1; I <= N; I++) {
        if (V > SEPTMP[I]) {
          TOL = ONE;
        } else {
          TOL = V / SEPTMP[I];
        }
        if (V > SEPIN[I]) {
          TOLIN = ONE;
        } else {
          TOLIN = V / SEPIN[I];
        }
        TOL = max(TOL, SMLNUM / EPS);
        TOLIN = max(TOLIN, SMLNUM / EPS);
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
        if (VMAX > RMAX[2]) {
          RMAX[2] = VMAX;
          if (NINFO[2] == 0) LMAX[2] = KNT.value;
        }
      }

      // Compare condition numbers for eigenvectors
      // taking their condition numbers into account

      for (I = 1; I <= N; I++) {
        if (V > SEPTMP[I] * STMP[I]) {
          TOL = SEPTMP[I];
        } else {
          TOL = V / STMP[I];
        }
        if (V > SEPIN[I] * SIN[I]) {
          TOLIN = SEPIN[I];
        } else {
          TOLIN = V / SIN[I];
        }
        TOL = max(TOL, SMLNUM / EPS);
        TOLIN = max(TOLIN, SMLNUM / EPS);
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
        if (VMAX > RMAX[2]) {
          RMAX[2] = VMAX;
          if (NINFO[2] == 0) LMAX[2] = KNT.value;
        }
      }

      // Compare condition numbers for eigenvalues
      // without taking their condition numbers into account

      for (I = 1; I <= N; I++) {
        if (SIN[I] <= 2 * N * EPS && STMP[I] <= 2 * N * EPS) {
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

      for (I = 1; I <= N; I++) {
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

      VMAX = ZERO;
      DUM[1] = -ONE;
      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      ztrsna('E', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != S[I]) VMAX = ONE / EPS;
        if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
      }

      // Compute eigenvector condition numbers only and compare

      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      ztrsna('V', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
        if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
      }

      // Compute all condition numbers using SELECT and compare

      for (I = 1; I <= N; I++) {
        SELECT[I] = true;
      }
      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      ztrsna('B', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
        if (STMP[I] != S[I]) VMAX = ONE / EPS;
      }

      // Compute eigenvalue condition numbers using SELECT and compare

      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      ztrsna('E', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != S[I]) VMAX = ONE / EPS;
        if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
      }

      // Compute eigenvector condition numbers using SELECT and compare

      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      ztrsna('V', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
        if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
      }
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }

      // Select second and next to last eigenvalues

      for (I = 1; I <= N; I++) {
        SELECT[I] = false;
      }
      ICMP = 0;
      if (N > 1) {
        ICMP = 1;
        LCMP[1] = 2;
        SELECT[2] = true;
        zcopy(N, RE(1, 2).asArray(), 1, RE(1, 1).asArray(), 1);
        zcopy(N, LE(1, 2).asArray(), 1, LE(1, 1).asArray(), 1);
      }
      if (N > 3) {
        ICMP = 2;
        LCMP[2] = N - 1;
        SELECT[N - 1] = true;
        zcopy(N, RE(1, N - 1).asArray(), 1, RE(1, 2).asArray(), 1);
        zcopy(N, LE(1, N - 1).asArray(), 1, LE(1, 2).asArray(), 1);
      }

      // Compute all selected condition numbers

      dcopy(ICMP, DUM, 0, STMP, 1);
      dcopy(ICMP, DUM, 0, SEPTMP, 1);
      ztrsna('B', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= ICMP; I++) {
        J = LCMP[I];
        if (SEPTMP[I] != SEP[J]) VMAX = ONE / EPS;
        if (STMP[I] != S[J]) VMAX = ONE / EPS;
      }

      // Compute selected eigenvalue condition numbers

      dcopy(ICMP, DUM, 0, STMP, 1);
      dcopy(ICMP, DUM, 0, SEPTMP, 1);
      ztrsna('E', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= ICMP; I++) {
        J = LCMP[I];
        if (STMP[I] != S[J]) VMAX = ONE / EPS;
        if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
      }

      // Compute selected eigenvector condition numbers

      dcopy(ICMP, DUM, 0, STMP, 1);
      dcopy(ICMP, DUM, 0, SEPTMP, 1);
      ztrsna('V', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M,
          WORK.asMatrix(), N, RWORK, INFO);
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3]++;
        continue;
      }
      for (I = 1; I <= ICMP; I++) {
        J = LCMP[I];
        if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
        if (SEPTMP[I] != SEP[J]) VMAX = ONE / EPS;
      }
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }
    }
  }
}
