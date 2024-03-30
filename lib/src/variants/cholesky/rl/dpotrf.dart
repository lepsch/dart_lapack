import 'dart:math';

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
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        final JB = min(NB, N - J + 1);

        dpotrf2('Upper', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }

        if (J + JB <= N) {
          // Updating the trailing submatrix.

          dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', JB, N - J - JB + 1,
              ONE, A(J, J), LDA, A(J, J + JB), LDA);
          dsyrk('Upper', 'Transpose', N - J - JB + 1, JB, -ONE, A(J, J + JB),
              LDA, ONE, A(J + JB, J + JB), LDA);
        }
      }
    } else {
      // Compute the Cholesky factorization A = L*L'.

      for (var J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        final JB = min(NB, N - J + 1);

        dpotrf2('Lower', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }

        if (J + JB <= N) {
          // Updating the trailing submatrix.

          dtrsm('Right', 'Lower', 'Transpose', 'Non-unit', N - J - JB + 1, JB,
              ONE, A(J, J), LDA, A(J + JB, J), LDA);
          dsyrk('Lower', 'No Transpose', N - J - JB + 1, JB, -ONE, A(J + JB, J),
              LDA, ONE, A(J + JB, J + JB), LDA);
        }
      }
    }
  }
}
