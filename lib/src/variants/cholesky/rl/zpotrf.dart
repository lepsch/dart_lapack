// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void zpotrf(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
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
    xerbla('ZPOTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine the block size for this environment.

  final NB = ilaenv(1, 'ZPOTRF', UPLO, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code.

    zpotrf2(UPLO, N, A, LDA, INFO);
  } else {
    // Use blocked code.

    if (UPPER) {
      // Compute the Cholesky factorization A = U'*U.

      for (var J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        final JB = min(NB, N - J + 1);

        zpotrf2('Upper', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }

        if (J + JB <= N) {
          // Updating the trailing submatrix.

          ztrsm('Left', 'Upper', 'Conjugate Transpose', 'Non-unit', JB,
              N - J - JB + 1, Complex.one, A(J, J), LDA, A(J, J + JB), LDA);
          zherk('Upper', 'Conjugate transpose', N - J - JB + 1, JB, -ONE,
              A(J, J + JB), LDA, ONE, A(J + JB, J + JB), LDA);
        }
      }
    } else {
      // Compute the Cholesky factorization A = L*L'.

      for (var J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        final JB = min(NB, N - J + 1);

        zpotrf2('Lower', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }

        if (J + JB <= N) {
          // Updating the trailing submatrix.

          ztrsm('Right', 'Lower', 'Conjugate Transpose', 'Non-unit',
              N - J - JB + 1, JB, Complex.one, A(J, J), LDA, A(J + JB, J), LDA);
          zherk('Lower', 'No Transpose', N - J - JB + 1, JB, -ONE, A(J + JB, J),
              LDA, ONE, A(J + JB, J + JB), LDA);
        }
      }
    }
  }
}
