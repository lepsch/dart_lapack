// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';

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
      // Compute the Cholesky factorization A = U'*U.

      for (var J = 1; J <= N; J += NB) {
        final JB = min(NB, N - J + 1);

        // Compute the current block.

        dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', J - 1, JB, ONE, A(1, 1),
            LDA, A(1, J), LDA);
        dsyrk('Upper', 'Transpose', JB, J - 1, -ONE, A(1, J), LDA, ONE, A(J, J),
            LDA);

        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        dpotrf2('Upper', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }
      }
    } else {
      // Compute the Cholesky factorization A = L*L'.

      for (var J = 1; J <= N; J += NB) {
        final JB = min(NB, N - J + 1);

        // Compute the current block.

        dtrsm('Right', 'Lower', 'Transpose', 'Non-unit', JB, J - 1, ONE,
            A(1, 1), LDA, A(J, 1), LDA);
        dsyrk('Lower', 'No Transpose', JB, J - 1, -ONE, A(J, 1), LDA, ONE,
            A(J, J), LDA);

        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        dpotrf2('Lower', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }
      }
    }
  }
}
