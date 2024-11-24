// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarf.dart';

void zunmr2(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool LEFT, NOTRAN;
  int I, I1, I2, I3, MI = 0, NI = 0, NQ;
  Complex AII, TAUI;

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');

  // NQ is the order of Q

  if (LEFT) {
    NQ = M;
  } else {
    NQ = N;
  }
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > NQ) {
    INFO.value = -5;
  } else if (LDA < max(1, K)) {
    INFO.value = -7;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZUNMR2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) return;

  if ((LEFT && !NOTRAN || !LEFT && NOTRAN)) {
    I1 = 1;
    I2 = K;
    I3 = 1;
  } else {
    I1 = K;
    I2 = 1;
    I3 = -1;
  }

  if (LEFT) {
    NI = N;
  } else {
    MI = M;
  }

  for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
    if (LEFT) {
      // H(i) or H(i)**H is applied to C(1:m-k+i,1:n)

      MI = M - K + I;
    } else {
      // H(i) or H(i)**H is applied to C(1:m,1:n-k+i)

      NI = N - K + I;
    }

    // Apply H(i) or H(i)**H

    if (NOTRAN) {
      TAUI = TAU[I].conjugate();
    } else {
      TAUI = TAU[I];
    }
    zlacgv(NQ - K + I - 1, A(I, 1).asArray(), LDA);
    AII = A[I][NQ - K + I];
    A[I][NQ - K + I] = Complex.one;
    zlarf(SIDE, MI, NI, A(I, 1).asArray(), LDA, TAUI, C, LDC, WORK);
    A[I][NQ - K + I] = AII;
    zlacgv(NQ - K + I - 1, A(I, 1).asArray(), LDA);
  }
}
