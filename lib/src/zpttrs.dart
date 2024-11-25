// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zptts2.dart';

void zpttrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final D = D_.having();
  final B = B_.having(ld: LDB);
  final E = E_.having();
  bool UPPER;
  int IUPLO, J, JB, NB;

  // Test the input arguments.

  INFO.value = 0;
  UPPER = (UPLO == 'U' || UPLO == 'u');
  if (!UPPER && !(UPLO == 'L' || UPLO == 'l')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZPTTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // Determine the number of right-hand sides to solve at a time.

  if (NRHS == 1) {
    NB = 1;
  } else {
    NB = max(1, ilaenv(1, 'ZPTTRS', UPLO, N, NRHS, -1, -1));
  }

  // Decode UPLO

  if (UPPER) {
    IUPLO = 1;
  } else {
    IUPLO = 0;
  }

  if (NB >= NRHS) {
    zptts2(IUPLO, N, NRHS, D, E, B, LDB);
  } else {
    for (J = 1; J <= NRHS; J += NB) {
      JB = min(NRHS - J + 1, NB);
      zptts2(IUPLO, N, JB, D, E, B(1, J), LDB);
    }
  }
}
