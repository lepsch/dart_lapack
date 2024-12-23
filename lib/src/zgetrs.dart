// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlaswp.dart';

void zgetrs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  bool NOTRAN;

  // Test the input parameters.

  INFO.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZGETRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (NOTRAN) {
    // Solve A * X = B.

    // Apply row interchanges to the right hand sides.

    zlaswp(NRHS, B, LDB, 1, N, IPIV, 1);

    // Solve L*X = B, overwriting B with X.

    ztrsm('Left', 'Lower', 'No transpose', 'Unit', N, NRHS, Complex.one, A, LDA,
        B, LDB);

    // Solve U*X = B, overwriting B with X.

    ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, Complex.one, A,
        LDA, B, LDB);
  } else {
    // Solve A**T * X = B  or A**H * X = B.

    // Solve U**T *X = B or U**H *X = B, overwriting B with X.

    ztrsm('Left', 'Upper', TRANS, 'Non-unit', N, NRHS, Complex.one, A, LDA, B,
        LDB);

    // Solve L**T *X = B, or L**H *X = B overwriting B with X.

    ztrsm('Left', 'Lower', TRANS, 'Unit', N, NRHS, Complex.one, A, LDA, B, LDB);

    // Apply row interchanges to the solution vectors.

    zlaswp(NRHS, B, LDB, 1, N, IPIV, -1);
  }
}
