// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dpttrf.dart';
import 'package:dart_lapack/src/dpttrs.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dptsv(
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DPTSV', -INFO.value);
    return;
  }

  // Compute the L*D*L**T (or U**T*D*U) factorization of A.

  dpttrf(N, D, E, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    dpttrs(N, NRHS, D, E, B, LDB, INFO);
  }
}
