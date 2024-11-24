// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/zdrot.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfgp.dart';
import 'package:lapack/src/zunbdb5.dart';

void zunbdb4(
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
  final Array<Complex> PHANTOM_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X11 = X11_.having(ld: LDX11);
  final X21 = X21_.having(ld: LDX21);
  final WORK = WORK_.having();
  final TAUP1 = TAUP1_.having();
  final TAUP2 = TAUP2_.having();
  final TAUQ1 = TAUQ1_.having();
  final PHANTOM = PHANTOM_.having();
  final THETA = THETA_.having();
  final PHI = PHI_.having();
  double C = 0, S = 0;
  int I, ILARF = 0, IORBDB5 = 0, J, LLARF, LORBDB5 = 0, LWORKMIN, LWORKOPT;
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
    WORK[1] = LWORKOPT.toComplex();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -14;
    }
  }
  if (INFO.value != 0) {
    xerbla('ZUNBDB4', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Reduce columns 1, ..., M-Q of X11 and X21

  for (I = 1; I <= M - Q; I++) {
    if (I == 1) {
      for (J = 1; J <= M; J++) {
        PHANTOM[J] = Complex.zero;
      }
      zunbdb5(P, M - P, Q, PHANTOM(1), 1, PHANTOM(P + 1), 1, X11, LDX11, X21,
          LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO);
      zscal(P, -Complex.one, PHANTOM(1), 1);
      zlarfgp(P, PHANTOM(1), PHANTOM(2), 1, TAUP1(1));
      zlarfgp(M - P, PHANTOM(P + 1), PHANTOM(P + 2), 1, TAUP2(1));
      THETA[I] = atan2(PHANTOM[1].real, PHANTOM[P + 1].real);
      C = cos(THETA[I]);
      S = sin(THETA[I]);
      PHANTOM[1] = Complex.one;
      PHANTOM[P + 1] = Complex.one;
      zlarf('L', P, Q, PHANTOM(1), 1, TAUP1[1].conjugate(), X11, LDX11,
          WORK(ILARF));
      zlarf('L', M - P, Q, PHANTOM(P + 1), 1, TAUP2[1].conjugate(), X21, LDX21,
          WORK(ILARF));
    } else {
      zunbdb5(
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
      zscal(P - I + 1, -Complex.one, X11(I, I - 1).asArray(), 1);
      zlarfgp(
          P - I + 1, X11(I, I - 1), X11(I + 1, I - 1).asArray(), 1, TAUP1(I));
      zlarfgp(M - P - I + 1, X21(I, I - 1), X21(I + 1, I - 1).asArray(), 1,
          TAUP2(I));
      THETA[I] = atan2(X11[I][I - 1].real, X21[I][I - 1].real);
      C = cos(THETA[I]);
      S = sin(THETA[I]);
      X11[I][I - 1] = Complex.one;
      X21[I][I - 1] = Complex.one;
      zlarf('L', P - I + 1, Q - I + 1, X11(I, I - 1).asArray(), 1,
          TAUP1[1].conjugate(), X11(I, I), LDX11, WORK(ILARF));
      zlarf('L', M - P - I + 1, Q - I + 1, X21(I, I - 1).asArray(), 1,
          TAUP2[1].conjugate(), X21(I, I), LDX21, WORK(ILARF));
    }

    zdrot(Q - I + 1, X11(I, I).asArray(), LDX11, X21(I, I).asArray(), LDX21, S,
        -C);
    zlacgv(Q - I + 1, X21(I, I).asArray(), LDX21);
    zlarfgp(Q - I + 1, X21(I, I), X21(I, I + 1).asArray(), LDX21, TAUQ1(I));
    C = X21[I][I].real;
    X21[I][I] = Complex.one;
    zlarf('R', P - I, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X11(I + 1, I), LDX11, WORK(ILARF));
    zlarf('R', M - P - I, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X21(I + 1, I), LDX21, WORK(ILARF));
    zlacgv(Q - I + 1, X21(I, I).asArray(), LDX21);
    if (I < M - Q) {
      S = sqrt(pow(dznrm2(P - I, X11(I + 1, I).asArray(), 1), 2) +
          pow(dznrm2(M - P - I, X21(I + 1, I).asArray(), 1), 2));
      PHI[I] = atan2(S, C);
    }
  }

  // Reduce the bottom-right portion of X11 to [ I 0 ]

  for (I = M - Q + 1; I <= P; I++) {
    zlacgv(Q - I + 1, X11(I, I).asArray(), LDX11);
    zlarfgp(Q - I + 1, X11(I, I), X11(I, I + 1).asArray(), LDX11, TAUQ1(I));
    X11[I][I] = Complex.one;
    zlarf('R', P - I, Q - I + 1, X11(I, I).asArray(), LDX11, TAUQ1[I],
        X11(I + 1, I), LDX11, WORK(ILARF));
    zlarf('R', Q - P, Q - I + 1, X11(I, I).asArray(), LDX11, TAUQ1[I],
        X21(M - Q + 1, I), LDX21, WORK(ILARF));
    zlacgv(Q - I + 1, X11(I, I).asArray(), LDX11);
  }

  // Reduce the bottom-right portion of X21 to [ 0 I ]

  for (I = P + 1; I <= Q; I++) {
    zlacgv(Q - I + 1, X21(M - Q + I - P, I).asArray(), LDX21);
    zlarfgp(Q - I + 1, X21(M - Q + I - P, I),
        X21(M - Q + I - P, I + 1).asArray(), LDX21, TAUQ1(I));
    X21[M - Q + I - P][I] = Complex.one;
    zlarf('R', Q - I, Q - I + 1, X21(M - Q + I - P, I).asArray(), LDX21,
        TAUQ1[I], X21(M - Q + I - P + 1, I), LDX21, WORK(ILARF));
    zlacgv(Q - I + 1, X21(M - Q + I - P, I).asArray(), LDX21);
  }
}
