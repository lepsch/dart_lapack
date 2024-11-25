// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void ztrti2(
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  bool NOUNIT, UPPER;
  int J;
  Complex AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOUNIT = lsame(DIAG, 'N');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZTRTI2', -INFO.value);
    return;
  }

  if (UPPER) {
    // Compute inverse of upper triangular matrix.

    for (J = 1; J <= N; J++) {
      if (NOUNIT) {
        A[J][J] = Complex.one / A[J][J];
        AJJ = -A[J][J];
      } else {
        AJJ = -Complex.one;
      }

      // Compute elements 1:j-1 of j-th column.

      ztrmv('Upper', 'No transpose', DIAG, J - 1, A, LDA, A(1, J).asArray(), 1);
      zscal(J - 1, AJJ, A(1, J).asArray(), 1);
    }
  } else {
    // Compute inverse of lower triangular matrix.

    for (J = N; J >= 1; J--) {
      if (NOUNIT) {
        A[J][J] = Complex.one / A[J][J];
        AJJ = -A[J][J];
      } else {
        AJJ = -Complex.one;
      }
      if (J < N) {
        // Compute elements j+1:n of j-th column.

        ztrmv('Lower', 'No transpose', DIAG, N - J, A(J + 1, J + 1), LDA,
            A(J + 1, J).asArray(), 1);
        zscal(N - J, AJJ, A(J + 1, J).asArray(), 1);
      }
    }
  }
}
