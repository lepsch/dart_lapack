import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgetf2.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgetrf(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  const ONE = 1.0;
  int I, J, JB, K, NB;
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
    xerbla('DGETRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'DGETRF', ' ', M, N, -1, -1);
  if (NB <= 1 || NB >= min(M, N)) {
    // Use unblocked code.

    dgetf2(M, N, A, LDA, IPIV, INFO);
  } else {
    // Use blocked code.

    for (J = 1; NB < 0 ? J >= min(M, N) : J <= min(M, N); J += NB) {
      JB = min(min(M, N) - J + 1, NB);

      // Update before factoring the current panel

      for (K = 1; NB < 0 ? K >= J - NB : K <= J - NB; K += NB) {
        // Apply interchanges to rows K:K+NB-1.

        dlaswp(JB, A(1, J), LDA, K, K + NB - 1, IPIV, 1);

        // Compute block row of U.

        dtrsm('Left', 'Lower', 'No transpose', 'Unit', NB, JB, ONE, A(K, K),
            LDA, A(K, J), LDA);

        // Update trailing submatrix.

        dgemm('No transpose', 'No transpose', M - K - NB + 1, JB, NB, -ONE,
            A(K + NB, K), LDA, A(K, J), LDA, ONE, A(K + NB, J), LDA);
      }

      // Factor diagonal and subdiagonal blocks and test for exact
      // singularity.

      dgetf2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO);

      // Adjust INFO.value and the pivot indices.

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + J - 1;
      for (I = J; I <= min(M, J + JB - 1); I++) {
        IPIV[I] = J - 1 + IPIV[I];
      }
    }

    // Apply interchanges to the left-overs

    for (K = 1; NB < 0 ? K >= min(M, N) : K <= min(M, N); K += NB) {
      dlaswp(K - 1, A(1, 1), LDA, K, min(K + NB - 1, min(M, N)), IPIV, 1);
    }

    // Apply update to the M+1:N columns when N > M

    if (N > M) {
      dlaswp(N - M, A(1, M + 1), LDA, 1, M, IPIV, 1);

      for (K = 1; NB < 0 ? K >= M : K <= M; K += NB) {
        JB = min(M - K + 1, NB);

        dtrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N - M, ONE, A(K, K),
            LDA, A(K, M + 1), LDA);

        if (K + NB <= M) {
          dgemm('No transpose', 'No transpose', M - K - NB + 1, N - M, NB, -ONE,
              A(K + NB, K), LDA, A(K, M + 1), LDA, ONE, A(K + NB, M + 1), LDA);
        }
      }
    }
  }
}
