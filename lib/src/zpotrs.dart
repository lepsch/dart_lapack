// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zpotrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  bool UPPER;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZPOTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B where A = U**H *U.

    // Solve U**H *X = B, overwriting B with X.

    ztrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', N, NRHS,
        Complex.one, A, LDA, B, LDB);

    // Solve U*X = B, overwriting B with X.

    ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, Complex.one, A,
        LDA, B, LDB);
  } else {
    // Solve A*X = B where A = L*L**H.

    // Solve L*X = B, overwriting B with X.

    ztrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, Complex.one, A,
        LDA, B, LDB);

    // Solve L**H *X = B, overwriting B with X.

    ztrsm('Left', 'Lower', 'Conjugate transpose', 'Non-unit', N, NRHS,
        Complex.one, A, LDA, B, LDB);
  }
}
