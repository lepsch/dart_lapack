import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/ztrti2.dart';

void ztrtri(
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
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
    xerbla('ZTRTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check for singularity if non-unit.

  if (NOUNIT) {
    for (INFO.value = 1; INFO.value <= N; INFO.value++) {
      // 10
      if (A[INFO.value][INFO.value] == Complex.zero) return;
    } // 10
    INFO.value = 0;
  }

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZTRTRI', UPLO + DIAG, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    ztrti2(UPLO, DIAG, N, A, LDA, INFO);
  } else {
    // Use blocked code

    if (UPPER) {
      // Compute inverse of upper triangular matrix

      for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) {
        // 20
        JB = min(NB, N - J + 1);

        // Compute rows 1:j-1 of current block column

        ztrmm('Left', 'Upper', 'No transpose', DIAG, J - 1, JB, Complex.one, A,
            LDA, A(1, J), LDA);
        ztrsm('Right', 'Upper', 'No transpose', DIAG, J - 1, JB, -Complex.one,
            A(J, J), LDA, A(1, J), LDA);

        // Compute inverse of current diagonal block

        ztrti2('Upper', DIAG, JB, A(J, J), LDA, INFO);
      } // 20
    } else {
      // Compute inverse of lower triangular matrix

      NN = ((N - 1) ~/ NB) * NB + 1;
      for (J = NN; -NB < 0 ? J >= 1 : J <= 1; J += -NB) {
        // 30
        JB = min(NB, N - J + 1);
        if (J + JB <= N) {
          // Compute rows j+jb:n of current block column

          ztrmm('Left', 'Lower', 'No transpose', DIAG, N - J - JB + 1, JB,
              Complex.one, A(J + JB, J + JB), LDA, A(J + JB, J), LDA);
          ztrsm('Right', 'Lower', 'No transpose', DIAG, N - J - JB + 1, JB,
              -Complex.one, A(J, J), LDA, A(J + JB, J), LDA);
        }

        // Compute inverse of current diagonal block

        ztrti2('Lower', DIAG, JB, A(J, J), LDA, INFO);
      } // 30
    }
  }
}
