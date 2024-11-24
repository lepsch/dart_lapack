// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void cgetrs(
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
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final IPIV = IPIV_.dim();
  final B = B_.dim(LDB);
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
    xerbla('CGETRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (NOTRAN) {
    // Solve A * X = B.

    // Apply row interchanges to the right hand sides.

    claswp(NRHS, B, LDB, 1, N, IPIV, 1);

    // Solve L*X = B, overwriting B with X.

    ctrsm('Left', 'Lower', 'No transpose', 'Unit', N, NRHS, Complex.one, A, LDA,
        B, LDB);

    // Solve U*X = B, overwriting B with X.

    ctrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, Complex.one, A,
        LDA, B, LDB);
  } else {
    // Solve A**T * X = B  or A**H * X = B.

    // Solve U**T *X = B or U**H *X = B, overwriting B with X.

    ctrsm('Left', 'Upper', TRANS, 'Non-unit', N, NRHS, Complex.one, A, LDA, B,
        LDB);

    // Solve L**T *X = B, or L**H *X = B overwriting B with X.

    ctrsm('Left', 'Lower', TRANS, 'Unit', N, NRHS, Complex.one, A, LDA, B, LDB);

    // Apply row interchanges to the solution vectors.

    claswp(NRHS, B, LDB, 1, N, IPIV, -1);
  }
}
