import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
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
  int J, JB, NB;

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
      // Compute the Cholesky factorization A = U**H *U.

      for (J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        JB = min(NB, N - J + 1);
        zherk('Upper', 'Conjugate transpose', JB, J - 1, -ONE, A(1, J), LDA,
            ONE, A(J, J), LDA);
        zpotrf2('Upper', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }
        if (J + JB <= N) {
          // Compute the current block row.

          zgemm(
              'Conjugate transpose',
              'No transpose',
              JB,
              N - J - JB + 1,
              J - 1,
              -Complex.one,
              A(1, J),
              LDA,
              A(1, J + JB),
              LDA,
              Complex.one,
              A(J, J + JB),
              LDA);
          ztrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', JB,
              N - J - JB + 1, Complex.one, A(J, J), LDA, A(J, J + JB), LDA);
        }
      }
    } else {
      // Compute the Cholesky factorization A = L*L**H.

      for (J = 1; J <= N; J += NB) {
        // Update and factorize the current diagonal block and test
        // for non-positive-definiteness.

        JB = min(NB, N - J + 1);
        zherk('Lower', 'No transpose', JB, J - 1, -ONE, A(J, 1), LDA, ONE,
            A(J, J), LDA);
        zpotrf2('Lower', JB, A(J, J), LDA, INFO);
        if (INFO.value != 0) {
          INFO.value += J - 1;
          return;
        }
        if (J + JB <= N) {
          // Compute the current block column.

          zgemm(
              'No transpose',
              'Conjugate transpose',
              N - J - JB + 1,
              JB,
              J - 1,
              -Complex.one,
              A(J + JB, 1),
              LDA,
              A(J, 1),
              LDA,
              Complex.one,
              A(J + JB, J),
              LDA);
          ztrsm('Right', 'Lower', 'Conjugate transpose', 'Non-unit',
              N - J - JB + 1, JB, Complex.one, A(J, J), LDA, A(J + JB, J), LDA);
        }
      }
    }
  }
}
