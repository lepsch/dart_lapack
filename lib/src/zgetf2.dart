// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/blas/zgeru.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zrscl.dart';

void zgetf2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  int J, JP;

  // Test the input parameters.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZGETF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  for (J = 1; J <= min(M, N); J++) {
    // Find pivot and test for singularity.

    JP = J - 1 + izamax(M - J + 1, A(J, J).asArray(), 1);
    IPIV[J] = JP;
    if (A[JP][J] != Complex.zero) {
      // Apply the interchange to columns 1:N.

      if (JP != J) zswap(N, A(J, 1).asArray(), LDA, A(JP, 1).asArray(), LDA);

      // Compute elements J+1:M of J-th column.

      if (J < M) zrscl(M - J, A[J][J], A(J + 1, J).asArray(), 1);
    } else if (INFO.value == 0) {
      INFO.value = J;
    }

    if (J < min(M, N)) {
      // Update trailing submatrix.

      zgeru(M - J, N - J, -Complex.one, A(J + 1, J).asArray(), 1,
          A(J, J + 1).asArray(), LDA, A(J + 1, J + 1), LDA);
    }
  }
}
