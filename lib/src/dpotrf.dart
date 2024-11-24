// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void dpotrf(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  const ONE = 1.0;

  // Test the input parameters.

  INFO.value = 0;
  final UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DPOTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine the block size for this environment.

  final NB = ilaenv(1, 'DPOTRF', UPLO, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code.

    dpotrf2(UPLO, N, A, LDA, INFO);
  } else {
    // Use blocked code.

    if (UPPER) {
      // Compute the Cholesky factorization A = U**T*U.

      for (var J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        final JB = min(NB, N - J + 1);
        dsyrk('Upper', 'Transpose', JB, J - 1, -ONE, A(1, J), LDA, ONE, A(J, J),
            LDA);
        dpotrf2('Upper', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }
        if (J + JB <= N) {
          // Compute the current block row.

          dgemm('Transpose', 'No transpose', JB, N - J - JB + 1, J - 1, -ONE,
              A(1, J), LDA, A(1, J + JB), LDA, ONE, A(J, J + JB), LDA);
          dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', JB, N - J - JB + 1,
              ONE, A(J, J), LDA, A(J, J + JB), LDA);
        }
      }
    } else {
      // Compute the Cholesky factorization A = L*L**T.

      for (var J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        final JB = min(NB, N - J + 1);
        dsyrk('Lower', 'No transpose', JB, J - 1, -ONE, A(J, 1), LDA, ONE,
            A(J, J), LDA);
        dpotrf2('Lower', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }
        if (J + JB <= N) {
          // Compute the current block column.

          dgemm('No transpose', 'Transpose', N - J - JB + 1, JB, J - 1, -ONE,
              A(J + JB, 1), LDA, A(J, 1), LDA, ONE, A(J + JB, J), LDA);
          dtrsm('Right', 'Lower', 'Transpose', 'Non-unit', N - J - JB + 1, JB,
              ONE, A(J, J), LDA, A(J + JB, J), LDA);
        }
      }
    }
  }
}
