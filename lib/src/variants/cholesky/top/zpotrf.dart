import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zpotrf2.dart';

void zpotrf(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  const ONE = 1.0;
  bool UPPER;
  int J = 0, JB, NB;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
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

  NB = ilaenv(1, 'ZPOTRF', UPLO, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code.

    zpotrf2(UPLO, N, A, LDA, INFO);
  } else {
    // Use blocked code.

    if (UPPER) {
      // Compute the Cholesky factorization A = U'*U.

      for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) {
        JB = min(NB, N - J + 1);

        // Compute the current block.

        ztrsm('Left', 'Upper', 'Conjugate Transpose', 'Non-unit', J - 1, JB,
            Complex.one, A(1, 1), LDA, A(1, J), LDA);
        zherk('Upper', 'Conjugate Transpose', JB, J - 1, -ONE, A(1, J), LDA,
            ONE, A(J, J), LDA);

        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        zpotrf2('Upper', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) break;
      }
    } else {
      // Compute the Cholesky factorization A = L*L'.

      for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) {
        JB = min(NB, N - J + 1);

        // Compute the current block.

        ztrsm('Right', 'Lower', 'Conjugate Transpose', 'Non-unit', JB, J - 1,
            Complex.one, A(1, 1), LDA, A(J, 1), LDA);
        zherk('Lower', 'No Transpose', JB, J - 1, -ONE, A(J, 1), LDA, ONE,
            A(J, J), LDA);

        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        zpotrf2('Lower', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) break;
      }
    }
  }

  if (INFO.value != 0) {
    INFO.value = INFO.value + J - 1;
  }
}
