// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlaorhr_col_getrfnp2.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlaorhr_col_getrfnp(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  const ONE = 1.0;
  int J, JB, NB;
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
    xerbla('DLAORHR_COL_GETRFNP', -INFO.value);
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'DLAORHR_COL_GETRFNP', ' ', M, N, -1, -1);

  if (NB <= 1 || NB >= min(M, N)) {
    // Use unblocked code.

    dlaorhr_col_getrfnp2(M, N, A, LDA, D, INFO);
  } else {
    // Use blocked code.

    for (J = 1; J <= min(M, N); J += NB) {
      JB = min(min(M, N) - J + 1, NB);

      // Factor diagonal and subdiagonal blocks.

      dlaorhr_col_getrfnp2(M - J + 1, JB, A(J, J), LDA, D(J), IINFO);

      if (J + JB <= N) {
        // Compute block row of U.

        dtrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N - J - JB + 1, ONE,
            A(J, J), LDA, A(J, J + JB), LDA);
        if (J + JB <= M) {
          // Update trailing submatrix.

          dgemm(
              'No transpose',
              'No transpose',
              M - J - JB + 1,
              N - J - JB + 1,
              JB,
              -ONE,
              A(J + JB, J),
              LDA,
              A(J, J + JB),
              LDA,
              ONE,
              A(J + JB, J + JB),
              LDA);
        }
      }
    }
  }
}
