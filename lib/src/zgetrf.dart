import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgetrf2.dart';
import 'package:lapack/src/zlaswp.dart';

void zgetrf(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  int I, J, JB, NB;
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
    xerbla('ZGETRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZGETRF', ' ', M, N, -1, -1);
  if (NB <= 1 || NB >= min(M, N)) {
    // Use unblocked code.

    zgetrf2(M, N, A, LDA, IPIV, INFO);
  } else {
    // Use blocked code.

    for (J = 1; NB < 0 ? J >= min(M, N) : J <= min(M, N); J += NB) {
      JB = min(min(M, N) - J + 1, NB);

      // Factor diagonal and subdiagonal blocks and test for exact
      // singularity.

      zgetrf2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO);

      // Adjust INFO.value and the pivot indices.

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + J - 1;
      for (I = J; I <= min(M, J + JB - 1); I++) {
        IPIV[I] = J - 1 + IPIV[I];
      }

      // Apply interchanges to columns 1:J-1.

      zlaswp(J - 1, A, LDA, J, J + JB - 1, IPIV, 1);

      if (J + JB <= N) {
        // Apply interchanges to columns J+JB:N.

        zlaswp(N - J - JB + 1, A(1, J + JB), LDA, J, J + JB - 1, IPIV, 1);

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
