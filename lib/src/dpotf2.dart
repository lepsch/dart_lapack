// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpotf2(
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
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int J;
  double AJJ;

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
    xerbla('DPOTF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Compute the Cholesky factorization A = U**T *U.

    for (J = 1; J <= N; J++) {
      // Compute U(J,J) and test for non-positive-definiteness.

      AJJ = A[J][J] - ddot(J - 1, A(1, J).asArray(), 1, A(1, J).asArray(), 1);
      if (AJJ <= ZERO || disnan(AJJ)) {
        A[J][J] = AJJ;
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      A[J][J] = AJJ;

      // Compute elements J+1:N of row J.

      if (J < N) {
        dgemv('Transpose', J - 1, N - J, -ONE, A(1, J + 1), LDA,
            A(1, J).asArray(), 1, ONE, A(J, J + 1).asArray(), LDA);
        dscal(N - J, ONE / AJJ, A(J, J + 1).asArray(), LDA);
      }
    }
  } else {
    // Compute the Cholesky factorization A = L*L**T.

    for (J = 1; J <= N; J++) {
      // Compute L(J,J) and test for non-positive-definiteness.

      AJJ =
          A[J][J] - ddot(J - 1, A(J, 1).asArray(), LDA, A(J, 1).asArray(), LDA);
      if (AJJ <= ZERO || disnan(AJJ)) {
        A[J][J] = AJJ;
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      A[J][J] = AJJ;

      // Compute elements J+1:N of column J.

      if (J < N) {
        dgemv('No transpose', N - J, J - 1, -ONE, A(J + 1, 1), LDA,
            A(J, 1).asArray(), LDA, ONE, A(J + 1, J).asArray(), 1);
        dscal(N - J, ONE / AJJ, A(J + 1, J).asArray(), 1);
      }
    }
  }
}
