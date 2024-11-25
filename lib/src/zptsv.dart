// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zpttrf.dart';
import 'package:dart_lapack/src/zpttrs.dart';

void zptsv(
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final B = B_.having(ld: LDB);
  final E = E_.having();
  final D = D_.having();

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
    xerbla('ZPTSV', -INFO.value);
    return;
  }

  // Compute the L*D*L**H (or U**H*D*U) factorization of A.

  zpttrf(N, D, E, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    zpttrs('Lower', N, NRHS, D, E, B, LDB, INFO);
  }
}
