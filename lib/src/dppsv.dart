// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dpptrf.dart';
import 'package:dart_lapack/src/dpptrs.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dppsv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> AP_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DPPSV', -INFO.value);
    return;
  }

  // Compute the Cholesky factorization A = U**T*U or A = L*L**T.

  dpptrf(UPLO, N, AP, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    dpptrs(UPLO, N, NRHS, AP, B, LDB, INFO);
  }
}
