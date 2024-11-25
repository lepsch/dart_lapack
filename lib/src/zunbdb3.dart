// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dznrm2.dart';
import 'package:dart_lapack/src/blas/zdrot.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';
import 'package:dart_lapack/src/zlarf.dart';
import 'package:dart_lapack/src/zlarfgp.dart';
import 'package:dart_lapack/src/zunbdb5.dart';

void zunbdb3(
  final int M,
  final int P,
  final int Q,
  final Matrix<Complex> X11_,
  final int LDX11,
  final Matrix<Complex> X21_,
  final int LDX21,
  final Array<double> THETA_,
  final Array<double> PHI_,
  final Array<Complex> TAUP1_,
  final Array<Complex> TAUP2_,
  final Array<Complex> TAUQ1_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final X11 = X11_.having(ld: LDX11);
  final X21 = X21_.having(ld: LDX21);
  final WORK = WORK_.having();
  final TAUP1 = TAUP1_.having();
  final TAUP2 = TAUP2_.having();
  final TAUQ1 = TAUQ1_.having();
  final THETA = THETA_.having();
  final PHI = PHI_.having();
  double C = 0, S = 0;
  int I, ILARF = 0, IORBDB5 = 0, LLARF, LORBDB5 = 0, LWORKMIN, LWORKOPT;
  bool LQUERY;
  final CHILDINFO = Box(0);

  // Test input arguments

  INFO.value = 0;
  LQUERY = LWORK == -1;

  if (M < 0) {
    INFO.value = -1;
  } else if (2 * P < M || P > M) {
    INFO.value = -2;
  } else if (Q < M - P || M - Q < M - P) {
    INFO.value = -3;
  } else if (LDX11 < max(1, P)) {
    INFO.value = -5;
  } else if (LDX21 < max(1, M - P)) {
    INFO.value = -7;
  }

  // Compute workspace

  if (INFO.value == 0) {
    ILARF = 2;
    LLARF = max(P, max(M - P - 1, Q - 1));
    IORBDB5 = 2;
    LORBDB5 = Q - 1;
    LWORKOPT = max(ILARF + LLARF - 1, IORBDB5 + LORBDB5 - 1);
    LWORKMIN = LWORKOPT;
    WORK[1] = LWORKOPT.toComplex();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -14;
    }
  }
  if (INFO.value != 0) {
    xerbla('ZUNBDB3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Reduce rows 1, ..., M-P of X11 and X21

  for (I = 1; I <= M - P; I++) {
    if (I > 1) {
      zdrot(Q - I + 1, X11(I - 1, I).asArray(), LDX11, X21(I, I).asArray(),
          LDX11, C, S);
    }

    zlacgv(Q - I + 1, X21(I, I).asArray(), LDX21);
    zlarfgp(Q - I + 1, X21(I, I), X21(I, I + 1).asArray(), LDX21, TAUQ1(I));
    S = X21[I][I].real;
    X21[I][I] = Complex.one;
    zlarf('R', P - I + 1, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X11(I, I), LDX11, WORK(ILARF));
    zlarf('R', M - P - I, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X21(I + 1, I), LDX21, WORK(ILARF));
    zlacgv(Q - I + 1, X21(I, I).asArray(), LDX21);
    C = sqrt(pow(dznrm2(P - I + 1, X11(I, I).asArray(), 1), 2) +
        pow(dznrm2(M - P - I, X21(I + 1, I).asArray(), 1), 2));
    THETA[I] = atan2(S, C);

    zunbdb5(
        P - I + 1,
        M - P - I,
        Q - I,
        X11(I, I).asArray(),
        1,
        X21(I + 1, I).asArray(),
        1,
        X11(I, I + 1),
        LDX11,
        X21(I + 1, I + 1),
        LDX21,
        WORK(IORBDB5),
        LORBDB5,
        CHILDINFO);
    zlarfgp(P - I + 1, X11(I, I), X11(I + 1, I).asArray(), 1, TAUP1(I));
    if (I < M - P) {
      zlarfgp(M - P - I, X21(I + 1, I), X21(I + 2, I).asArray(), 1, TAUP2(I));
      PHI[I] = atan2(X21[I + 1][I].real, X11[I][I].real);
      C = cos(PHI[I]);
      S = sin(PHI[I]);
      X21[I + 1][I] = Complex.one;
      zlarf('L', M - P - I, Q - I, X21(I + 1, I).asArray(), 1,
          TAUP2[I].conjugate(), X21(I + 1, I + 1), LDX21, WORK(ILARF));
    }
    X11[I][I] = Complex.one;
    zlarf('L', P - I + 1, Q - I, X11(I, I).asArray(), 1, TAUP1[I].conjugate(),
        X11(I, I + 1), LDX11, WORK(ILARF));
  }

  // Reduce the bottom-right portion of X11 to the identity matrix

  for (I = M - P + 1; I <= Q; I++) {
    zlarfgp(P - I + 1, X11(I, I), X11(I + 1, I).asArray(), 1, TAUP1(I));
    X11[I][I] = Complex.one;
    zlarf('L', P - I + 1, Q - I, X11(I, I).asArray(), 1, TAUP1[I].conjugate(),
        X11(I, I + 1), LDX11, WORK(ILARF));
  }
}
