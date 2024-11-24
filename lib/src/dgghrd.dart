// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgghrd(
  final String COMPQ,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  const ONE = 1.0, ZERO = 0.0;
  bool ILQ = false, ILZ = false;
  int ICOMPQ, ICOMPZ, JCOL, JROW;
  double TEMP;
  final C = Box(0.0), S = Box(0.0);

  // Decode COMPQ

  if (lsame(COMPQ, 'N')) {
    ILQ = false;
    ICOMPQ = 1;
  } else if (lsame(COMPQ, 'V')) {
    ILQ = true;
    ICOMPQ = 2;
  } else if (lsame(COMPQ, 'I')) {
    ILQ = true;
    ICOMPQ = 3;
  } else {
    ICOMPQ = 0;
  }

  // Decode COMPZ

  if (lsame(COMPZ, 'N')) {
    ILZ = false;
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'V')) {
    ILZ = true;
    ICOMPZ = 2;
  } else if (lsame(COMPZ, 'I')) {
    ILZ = true;
    ICOMPZ = 3;
  } else {
    ICOMPZ = 0;
  }

  // Test the input parameters.

  INFO.value = 0;
  if (ICOMPQ <= 0) {
    INFO.value = -1;
  } else if (ICOMPZ <= 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ILO < 1) {
    INFO.value = -4;
  } else if (IHI > N || IHI < ILO - 1) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if ((ILQ && LDQ < N) || LDQ < 1) {
    INFO.value = -11;
  } else if ((ILZ && LDZ < N) || LDZ < 1) {
    INFO.value = -13;
  }
  if (INFO.value != 0) {
    xerbla('DGGHRD', -INFO.value);
    return;
  }

  // Initialize Q and Z if desired.

  if (ICOMPQ == 3) dlaset('Full', N, N, ZERO, ONE, Q, LDQ);
  if (ICOMPZ == 3) dlaset('Full', N, N, ZERO, ONE, Z, LDZ);

  // Quick return if possible

  if (N <= 1) return;

  // Zero out lower triangle of B

  for (JCOL = 1; JCOL <= N - 1; JCOL++) {
    for (JROW = JCOL + 1; JROW <= N; JROW++) {
      B[JROW][JCOL] = ZERO;
    }
  }

  // Reduce A and B

  for (JCOL = ILO; JCOL <= IHI - 2; JCOL++) {
    for (JROW = IHI; JROW >= JCOL + 2; JROW--) {
      // Step 1: rotate rows JROW-1, JROW to kill A[JROW][JCOL]

      TEMP = A[JROW - 1][JCOL];
      dlartg(TEMP, A[JROW][JCOL], C, S, A.box(JROW - 1, JCOL));
      A[JROW][JCOL] = ZERO;
      drot(N - JCOL, A(JROW - 1, JCOL + 1).asArray(), LDA,
          A(JROW, JCOL + 1).asArray(), LDA, C.value, S.value);
      drot(N + 2 - JROW, B(JROW - 1, JROW - 1).asArray(), LDB,
          B(JROW, JROW - 1).asArray(), LDB, C.value, S.value);
      if (ILQ) {
        drot(N, Q(1, JROW - 1).asArray(), 1, Q(1, JROW).asArray(), 1, C.value,
            S.value);
      }

      // Step 2: rotate columns JROW, JROW-1 to kill B[JROW][JROW-1]

      TEMP = B[JROW][JROW];
      dlartg(TEMP, B[JROW][JROW - 1], C, S, B.box(JROW, JROW));
      B[JROW][JROW - 1] = ZERO;
      drot(IHI, A(1, JROW).asArray(), 1, A(1, JROW - 1).asArray(), 1, C.value,
          S.value);
      drot(JROW - 1, B(1, JROW).asArray(), 1, B(1, JROW - 1).asArray(), 1,
          C.value, S.value);
      if (ILZ) {
        drot(N, Z(1, JROW).asArray(), 1, Z(1, JROW - 1).asArray(), 1, C.value,
            S.value);
      }
    }
  }
}
