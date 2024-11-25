// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarfgp.dart';
import 'package:dart_lapack/src/dorbdb5.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorbdb4(
  final int M,
  final int P,
  final int Q,
  final Matrix<double> X11_,
  final int LDX11,
  final Matrix<double> X21_,
  final int LDX21,
  final Array<double> THETA_,
  final Array<double> PHI_,
  final Array<double> TAUP1_,
  final Array<double> TAUP2_,
  final Array<double> TAUQ1_,
  final Array<double> PHANTOM_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final X11 = X11_.having(ld: LDX11);
  final X21 = X21_.having(ld: LDX21);
  final THETA = THETA_.having();
  final PHI = PHI_.having();
  final TAUP1 = TAUP1_.having();
  final TAUP2 = TAUP2_.having();
  final TAUQ1 = TAUQ1_.having();
  final PHANTOM = PHANTOM_.having();
  final WORK = WORK_.having();
  const NEGONE = -1.0, ONE = 1.0, ZERO = 0.0;
  double C, S;
  int I, ILARF = 0, IORBDB5 = 0, J, LLARF, LORBDB5 = 0, LWORKMIN = 0, LWORKOPT;
  bool LQUERY;
  final CHILDINFO = Box(0);

  // Test input arguments

  INFO.value = 0;
  LQUERY = LWORK == -1;

  if (M < 0) {
    INFO.value = -1;
  } else if (P < M - Q || M - P < M - Q) {
    INFO.value = -2;
  } else if (Q < M - Q || Q > M) {
    INFO.value = -3;
  } else if (LDX11 < max(1, P)) {
    INFO.value = -5;
  } else if (LDX21 < max(1, M - P)) {
    INFO.value = -7;
  }

  // Compute workspace

  if (INFO.value == 0) {
    ILARF = 2;
    LLARF = max(Q - 1, max(P - 1, M - P - 1));
    IORBDB5 = 2;
    LORBDB5 = Q;
    LWORKOPT = ILARF + LLARF - 1;
    LWORKOPT = max(LWORKOPT, IORBDB5 + LORBDB5 - 1);
    LWORKMIN = LWORKOPT;
    WORK[1] = LWORKOPT.toDouble();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -14;
    }
  }
  if (INFO.value != 0) {
    xerbla('DORBDB4', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Reduce columns 1, ..., M-Q of X11 and X21

  for (I = 1; I <= M - Q; I++) {
    if (I == 1) {
      for (J = 1; J <= M; J++) {
        PHANTOM[J] = ZERO;
      }
      dorbdb5(P, M - P, Q, PHANTOM(1), 1, PHANTOM(P + 1), 1, X11, LDX11, X21,
          LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO);
      dscal(P, NEGONE, PHANTOM(1), 1);
      dlarfgp(P, PHANTOM.box(1), PHANTOM(2), 1, TAUP1.box(1));
      dlarfgp(M - P, PHANTOM.box(P + 1), PHANTOM(P + 2), 1, TAUP2.box(1));
      THETA[I] = atan2(PHANTOM[1], PHANTOM[P + 1]);
      C = cos(THETA[I]);
      S = sin(THETA[I]);
      PHANTOM[1] = ONE;
      PHANTOM[P + 1] = ONE;
      dlarf('L', P, Q, PHANTOM(1), 1, TAUP1[1], X11, LDX11, WORK(ILARF));
      dlarf(
          'L', M - P, Q, PHANTOM(P + 1), 1, TAUP2[1], X21, LDX21, WORK(ILARF));
    } else {
      dorbdb5(
          P - I + 1,
          M - P - I + 1,
          Q - I + 1,
          X11(I, I - 1).asArray(),
          1,
          X21(I, I - 1).asArray(),
          1,
          X11(I, I),
          LDX11,
          X21(I, I),
          LDX21,
          WORK(IORBDB5),
          LORBDB5,
          CHILDINFO);
      dscal(P - I + 1, NEGONE, X11(I, I - 1).asArray(), 1);
      dlarfgp(P - I + 1, X11.box(I, I - 1), X11(I + 1, I - 1).asArray(), 1,
          TAUP1.box(I));
      dlarfgp(M - P - I + 1, X21.box(I, I - 1), X21(I + 1, I - 1).asArray(), 1,
          TAUP2.box(I));
      THETA[I] = atan2(X11[I][I - 1], X21[I][I - 1]);
      C = cos(THETA[I]);
      S = sin(THETA[I]);
      X11[I][I - 1] = ONE;
      X21[I][I - 1] = ONE;
      dlarf('L', P - I + 1, Q - I + 1, X11(I, I - 1).asArray(), 1, TAUP1[I],
          X11(I, I), LDX11, WORK(ILARF));
      dlarf('L', M - P - I + 1, Q - I + 1, X21(I, I - 1).asArray(), 1, TAUP2[I],
          X21(I, I), LDX21, WORK(ILARF));
    }

    drot(Q - I + 1, X11(I, I).asArray(), LDX11, X21(I, I).asArray(), LDX21, S,
        -C);
    dlarfgp(
        Q - I + 1, X21.box(I, I), X21(I, I + 1).asArray(), LDX21, TAUQ1.box(I));
    C = X21[I][I];
    X21[I][I] = ONE;
    dlarf('R', P - I, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X11(I + 1, I), LDX11, WORK(ILARF));
    dlarf('R', M - P - I, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X21(I + 1, I), LDX21, WORK(ILARF));
    if (I < M - Q) {
      S = sqrt(pow(dnrm2(P - I, X11(I + 1, I).asArray(), 1), 2) +
          pow(dnrm2(M - P - I, X21(I + 1, I).asArray(), 1), 2));
      PHI[I] = atan2(S, C);
    }
  }

  // Reduce the bottom-right portion of X11 to [ I 0 ]

  for (I = M - Q + 1; I <= P; I++) {
    dlarfgp(
        Q - I + 1, X11.box(I, I), X11(I, I + 1).asArray(), LDX11, TAUQ1.box(I));
    X11[I][I] = ONE;
    dlarf('R', P - I, Q - I + 1, X11(I, I).asArray(), LDX11, TAUQ1[I],
        X11(I + 1, I), LDX11, WORK(ILARF));
    dlarf('R', Q - P, Q - I + 1, X11(I, I).asArray(), LDX11, TAUQ1[I],
        X21(M - Q + 1, I), LDX21, WORK(ILARF));
  }

  // Reduce the bottom-right portion of X21 to [ 0 I ]

  for (I = P + 1; I <= Q; I++) {
    dlarfgp(Q - I + 1, X21.box(M - Q + I - P, I),
        X21(M - Q + I - P, I + 1).asArray(), LDX21, TAUQ1.box(I));
    X21[M - Q + I - P][I] = ONE;
    dlarf('R', Q - I, Q - I + 1, X21(M - Q + I - P, I).asArray(), LDX21,
        TAUQ1[I], X21(M - Q + I - P + 1, I), LDX21, WORK(ILARF));
  }
}
