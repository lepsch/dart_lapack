// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgetf2.dart';
import 'package:dart_lapack/src/zlaswp.dart';

void zgetrf(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  int I, J, JB, NB;
  final IINFO = Box(0);

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
    xerbla('ZGETRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZGETRF', ' ', M, N, -1, -1);
  if (NB <= 1 || NB >= min(M, N)) {
    // Use unblocked code.

    zgetf2(M, N, A, LDA, IPIV, INFO);
  } else {
    // Use blocked code.

    for (J = 1; J <= min(M, N); J += NB) {
      JB = min(min(M, N) - J + 1, NB);

      // Update current block.

      zgemm('No transpose', 'No transpose', M - J + 1, JB, J - 1, -Complex.one,
          A(J, 1), LDA, A(1, J), LDA, Complex.one, A(J, J), LDA);

      // Factor diagonal and subdiagonal blocks and test for exact
      // singularity.

      zgetf2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO);

      // Adjust INFO and the pivot indices.

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + J - 1;
      for (I = J; I <= min(M, J + JB - 1); I++) {
        IPIV[I] = J - 1 + IPIV[I];
      }

      // Apply interchanges to column 1:J-1

      zlaswp(J - 1, A, LDA, J, J + JB - 1, IPIV, 1);

      if (J + JB <= N) {
        // Apply interchanges to column J+JB:N

        zlaswp(N - J - JB + 1, A(1, J + JB), LDA, J, J + JB - 1, IPIV, 1);

        zgemm(
            'No transpose',
            'No transpose',
            JB,
            N - J - JB + 1,
            J - 1,
            -Complex.one,
            A(J, 1),
            LDA,
            A(1, J + JB),
            LDA,
            Complex.one,
            A(J, J + JB),
            LDA);

        // Compute block row of U.

        ztrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N - J - JB + 1,
            Complex.one, A(J, J), LDA, A(J, J + JB), LDA);
      }
    }
  }
}
