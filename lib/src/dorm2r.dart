// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorm2r(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LEFT, NOTRAN;
  int I, I1, I2, I3, IC = 0, JC = 0, MI = 0, NI = 0, NQ;
  double AII;

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
  } else if (!NOTRAN && !lsame(TRANS, 'T')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > NQ) {
    INFO.value = -5;
  } else if (LDA < max(1, NQ)) {
    INFO.value = -7;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DORM2R', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) return;

  if ((LEFT && !NOTRAN) || (!LEFT && NOTRAN)) {
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
    JC = 1;
  } else {
    MI = M;
    IC = 1;
  }

  for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
    if (LEFT) {
      // H(i) is applied to C[i:m][1:n]

      MI = M - I + 1;
      IC = I;
    } else {
      // H(i) is applied to C[1:m][i:n]

      NI = N - I + 1;
      JC = I;
    }

    // Apply H(i)

    AII = A[I][I];
    A[I][I] = ONE;
    dlarf(SIDE, MI, NI, A(I, I).asArray(), 1, TAU[I], C(IC, JC), LDC, WORK);
    A[I][I] = AII;
  }
}
