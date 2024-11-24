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
import 'package:dart_lapack/src/zlaunhr_col_getrfnp2.dart';

void zlaunhr_col_getrfnp(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> D_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
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
    xerbla('ZLAUNHR_COL_GETRFNP', -INFO.value);
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZLAUNHR_COL_GETRFNP', ' ', M, N, -1, -1);

  if (NB <= 1 || NB >= min(M, N)) {
    // Use unblocked code.

    zlaunhr_col_getrfnp2(M, N, A, LDA, D, INFO);
  } else {
    // Use blocked code.

    for (J = 1; J <= min(M, N); J += NB) {
      JB = min(min(M, N) - J + 1, NB);

      // Factor diagonal and subdiagonal blocks.

      zlaunhr_col_getrfnp2(M - J + 1, JB, A(J, J), LDA, D(J), IINFO);

      if (J + JB <= N) {
        // Compute block row of U.

        ztrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N - J - JB + 1,
            Complex.one, A(J, J), LDA, A(J, J + JB), LDA);
        if (J + JB <= M) {
          // Update trailing submatrix.

          zgemm(
              'No transpose',
              'No transpose',
              M - J - JB + 1,
              N - J - JB + 1,
              JB,
              -Complex.one,
              A(J + JB, J),
              LDA,
              A(J, J + JB),
              LDA,
              Complex.one,
              A(J + JB, J + JB),
              LDA);
        }
      }
    }
  }
}
