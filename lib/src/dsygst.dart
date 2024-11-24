// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dsymm.dart';
import 'package:lapack/src/blas/dsyr2k.dart';
import 'package:lapack/src/blas/dtrmm.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dsygs2.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsygst(
  final int ITYPE,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, HALF = 0.5;
  bool UPPER;
  int K, KB, NB;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DSYGST', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'DSYGST', UPLO, N, -1, -1, -1);

  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    dsygs2(ITYPE, UPLO, N, A, LDA, B, LDB, INFO);
  } else {
    // Use blocked code

    if (ITYPE == 1) {
      if (UPPER) {
        // Compute inv(U**T)*A*inv(U)

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the upper triangle of A(k:n,k:n)

          dsygs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
          if (K + KB <= N) {
            dtrsm('Left', UPLO, 'Transpose', 'Non-unit', KB, N - K - KB + 1,
                ONE, B(K, K), LDB, A(K, K + KB), LDA);
            dsymm('Left', UPLO, KB, N - K - KB + 1, -HALF, A(K, K), LDA,
                B(K, K + KB), LDB, ONE, A(K, K + KB), LDA);
            dsyr2k(UPLO, 'Transpose', N - K - KB + 1, KB, -ONE, A(K, K + KB),
                LDA, B(K, K + KB), LDB, ONE, A(K + KB, K + KB), LDA);
            dsymm('Left', UPLO, KB, N - K - KB + 1, -HALF, A(K, K), LDA,
                B(K, K + KB), LDB, ONE, A(K, K + KB), LDA);
            dtrsm('Right', UPLO, 'No transpose', 'Non-unit', KB, N - K - KB + 1,
                ONE, B(K + KB, K + KB), LDB, A(K, K + KB), LDA);
          }
        }
      } else {
        // Compute inv(L)*A*inv(L**T)

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the lower triangle of A(k:n,k:n)

          dsygs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
          if (K + KB <= N) {
            dtrsm('Right', UPLO, 'Transpose', 'Non-unit', N - K - KB + 1, KB,
                ONE, B(K, K), LDB, A(K + KB, K), LDA);
            dsymm('Right', UPLO, N - K - KB + 1, KB, -HALF, A(K, K), LDA,
                B(K + KB, K), LDB, ONE, A(K + KB, K), LDA);
            dsyr2k(UPLO, 'No transpose', N - K - KB + 1, KB, -ONE, A(K + KB, K),
                LDA, B(K + KB, K), LDB, ONE, A(K + KB, K + KB), LDA);
            dsymm('Right', UPLO, N - K - KB + 1, KB, -HALF, A(K, K), LDA,
                B(K + KB, K), LDB, ONE, A(K + KB, K), LDA);
            dtrsm('Left', UPLO, 'No transpose', 'Non-unit', N - K - KB + 1, KB,
                ONE, B(K + KB, K + KB), LDB, A(K + KB, K), LDA);
          }
        }
      }
    } else {
      if (UPPER) {
        // Compute U*A*U**T

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the upper triangle of A(1:k+kb-1,1:k+kb-1)

          dtrmm('Left', UPLO, 'No transpose', 'Non-unit', K - 1, KB, ONE, B,
              LDB, A(1, K), LDA);
          dsymm('Right', UPLO, K - 1, KB, HALF, A(K, K), LDA, B(1, K), LDB, ONE,
              A(1, K), LDA);
          dsyr2k(UPLO, 'No transpose', K - 1, KB, ONE, A(1, K), LDA, B(1, K),
              LDB, ONE, A, LDA);
          dsymm('Right', UPLO, K - 1, KB, HALF, A(K, K), LDA, B(1, K), LDB, ONE,
              A(1, K), LDA);
          dtrmm('Right', UPLO, 'Transpose', 'Non-unit', K - 1, KB, ONE, B(K, K),
              LDB, A(1, K), LDA);
          dsygs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
        }
      } else {
        // Compute L**T*A*L

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the lower triangle of A(1:k+kb-1,1:k+kb-1)

          dtrmm('Right', UPLO, 'No transpose', 'Non-unit', KB, K - 1, ONE, B,
              LDB, A(K, 1), LDA);
          dsymm('Left', UPLO, KB, K - 1, HALF, A(K, K), LDA, B(K, 1), LDB, ONE,
              A(K, 1), LDA);
          dsyr2k(UPLO, 'Transpose', K - 1, KB, ONE, A(K, 1), LDA, B(K, 1), LDB,
              ONE, A, LDA);
          dsymm('Left', UPLO, KB, K - 1, HALF, A(K, K), LDA, B(K, 1), LDB, ONE,
              A(K, 1), LDA);
          dtrmm('Left', UPLO, 'Transpose', 'Non-unit', KB, K - 1, ONE, B(K, K),
              LDB, A(K, 1), LDA);
          dsygs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
        }
      }
    }
  }
}
