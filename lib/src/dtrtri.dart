// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtrmm.dart';
import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dtrti2.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtrtri(
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  const ONE = 1.0, ZERO = 0.0;
  bool NOUNIT, UPPER;
  int J, JB, NB, NN;

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
    xerbla('DTRTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check for singularity if non-unit.

  if (NOUNIT) {
    for (INFO.value = 1; INFO.value <= N; INFO.value++) {
      if (A[INFO.value][INFO.value] == ZERO) return;
    }
    INFO.value = 0;
  }

  // Determine the block size for this environment.

  NB = ilaenv(1, 'DTRTRI', UPLO + DIAG, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    dtrti2(UPLO, DIAG, N, A, LDA, INFO);
  } else {
    // Use blocked code

    if (UPPER) {
      // Compute inverse of upper triangular matrix

      for (J = 1; J <= N; J += NB) {
        JB = min(NB, N - J + 1);

        // Compute rows 1:j-1 of current block column

        dtrmm('Left', 'Upper', 'No transpose', DIAG, J - 1, JB, ONE, A, LDA,
            A(1, J), LDA);
        dtrsm('Right', 'Upper', 'No transpose', DIAG, J - 1, JB, -ONE, A(J, J),
            LDA, A(1, J), LDA);

        // Compute inverse of current diagonal block

        dtrti2('Upper', DIAG, JB, A(J, J), LDA, INFO);
      }
    } else {
      // Compute inverse of lower triangular matrix

      NN = ((N - 1) ~/ NB) * NB + 1;
      for (J = NN; J >= 1; J -= NB) {
        JB = min(NB, N - J + 1);
        if (J + JB <= N) {
          // Compute rows j+jb:n of current block column

          dtrmm('Left', 'Lower', 'No transpose', DIAG, N - J - JB + 1, JB, ONE,
              A(J + JB, J + JB), LDA, A(J + JB, J), LDA);
          dtrsm('Right', 'Lower', 'No transpose', DIAG, N - J - JB + 1, JB,
              -ONE, A(J, J), LDA, A(J + JB, J), LDA);
        }

        // Compute inverse of current diagonal block

        dtrti2('Lower', DIAG, JB, A(J, J), LDA, INFO);
      }
    }
  }
}
