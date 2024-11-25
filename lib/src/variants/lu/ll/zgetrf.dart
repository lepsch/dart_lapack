// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgetf2.dart';
import 'package:dart_lapack/src/zlaswp.dart';

void zgetrf(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
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
    xerbla('ZGETRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZGETRF', ' ', M, N, -1, -1);
  if (NB <= 1 || NB >= min(M, N)) {
    // Use unblocked code.

    zgetf2(M, N, A, LDA, IPIV, INFO);
  } else {
    // Use blocked code.

    for (J = 1; J <= min(M, N); J += NB) {
      JB = min(min(M, N) - J + 1, NB);

      // Update before factoring the current panel

      for (K = 1; K <= J - NB; K += NB) {
        // Apply interchanges to rows K:K+NB-1.

        zlaswp(JB, A(1, J), LDA, K, K + NB - 1, IPIV, 1);

        // Compute block row of U.

        ztrsm('Left', 'Lower', 'No transpose', 'Unit', NB, JB, Complex.one,
            A(K, K), LDA, A(K, J), LDA);

        // Update trailing submatrix.

        zgemm(
            'No transpose',
            'No transpose',
            M - K - NB + 1,
            JB,
            NB,
            -Complex.one,
            A(K + NB, K),
            LDA,
            A(K, J),
            LDA,
            Complex.one,
            A(K + NB, J),
            LDA);
      }

      // Factor diagonal and subdiagonal blocks and test for exact
      // singularity.

      zgetf2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO);

      // Adjust INFO and the pivot indices.

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + J - 1;
      for (I = J; I <= min(M, J + JB - 1); I++) {
        IPIV[I] = J - 1 + IPIV[I];
      }
    }

    // Apply interchanges to the left-overs

    for (K = 1; K <= min(M, N); K += NB) {
      zlaswp(K - 1, A(1, 1), LDA, K, min(K + NB - 1, min(M, N)), IPIV, 1);
    }

    // Apply update to the M+1:N columns when N > M

    if (N > M) {
      zlaswp(N - M, A(1, M + 1), LDA, 1, M, IPIV, 1);

      for (K = 1; K <= M; K += NB) {
        JB = min(M - K + 1, NB);

        ztrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N - M, Complex.one,
            A(K, K), LDA, A(K, M + 1), LDA);

        if (K + NB <= M) {
          zgemm(
              'No transpose',
              'No transpose',
              M - K - NB + 1,
              N - M,
              NB,
              -Complex.one,
              A(K + NB, K),
              LDA,
              A(K, M + 1),
              LDA,
              Complex.one,
              A(K + NB, M + 1),
              LDA);
        }
      }
    }
  }
}
