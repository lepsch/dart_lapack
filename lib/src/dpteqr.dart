// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dbdsqr.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dpttrf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dpteqr(
  final String COMPZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  final C = Matrix<double>(1, 1), VT = Matrix<double>(1, 1);
  int I, ICOMPZ, NRU;

  // Test the input parameters.

  INFO.value = 0;

  if (lsame(COMPZ, 'N')) {
    ICOMPZ = 0;
  } else if (lsame(COMPZ, 'V')) {
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'I')) {
    ICOMPZ = 2;
  } else {
    ICOMPZ = -1;
  }
  if (ICOMPZ < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if ((LDZ < 1) || (ICOMPZ > 0 && LDZ < max(1, N))) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DPTEQR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    if (ICOMPZ > 0) Z[1][1] = ONE;
    return;
  }
  if (ICOMPZ == 2) dlaset('Full', N, N, ZERO, ONE, Z, LDZ);

  // Call DPTTRF to factor the matrix.

  dpttrf(N, D, E, INFO);
  if (INFO.value != 0) return;
  for (I = 1; I <= N; I++) {
    D[I] = sqrt(D[I]);
  }
  for (I = 1; I <= N - 1; I++) {
    E[I] *= D[I];
  }

  // Call DBDSQR to compute the singular values/vectors of the
  // bidiagonal factor.

  if (ICOMPZ > 0) {
    NRU = N;
  } else {
    NRU = 0;
  }
  dbdsqr('Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1, WORK, INFO);

  // Square the singular values.

  if (INFO.value == 0) {
    for (I = 1; I <= N; I++) {
      D[I] *= D[I];
    }
  } else {
    INFO.value = N + INFO.value;
  }
}
