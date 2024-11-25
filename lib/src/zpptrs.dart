// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/ztpsv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zpptrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final B = B_.having(ld: LDB);
  bool UPPER;
  int I;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZPPTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B where A = U**H * U.

    for (I = 1; I <= NRHS; I++) {
      // Solve U**H *X = B, overwriting B with X.

      ztpsv('Upper', 'Conjugate transpose', 'Non-unit', N, AP,
          B(1, I).asArray(), 1);

      // Solve U*X = B, overwriting B with X.

      ztpsv('Upper', 'No transpose', 'Non-unit', N, AP, B(1, I).asArray(), 1);
    }
  } else {
    // Solve A*X = B where A = L * L**H.

    for (I = 1; I <= NRHS; I++) {
      // Solve L*Y = B, overwriting B with X.

      ztpsv('Lower', 'No transpose', 'Non-unit', N, AP, B(1, I).asArray(), 1);

      // Solve L**H *X = Y, overwriting B with X.

      ztpsv('Lower', 'Conjugate transpose', 'Non-unit', N, AP,
          B(1, I).asArray(), 1);
    }
  }
}
