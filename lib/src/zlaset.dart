// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlaset(
  final String UPLO,
  final int M,
  final int N,
  final Complex ALPHA,
  final Complex BETA,
  final Matrix<Complex> A_,
  final int LDA,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  int I, J;

  if (lsame(UPLO, 'U')) {
    // Set the diagonal to BETA and the strictly upper triangular
    // part of the array to ALPHA.

    for (J = 2; J <= N; J++) {
      for (I = 1; I <= min(J - 1, M); I++) {
        A[I][J] = ALPHA;
      }
    }
    for (I = 1; I <= min(N, M); I++) {
      A[I][I] = BETA;
    }
  } else if (lsame(UPLO, 'L')) {
    // Set the diagonal to BETA and the strictly lower triangular
    // part of the array to ALPHA.

    for (J = 1; J <= min(M, N); J++) {
      for (I = J + 1; I <= M; I++) {
        A[I][J] = ALPHA;
      }
    }
    for (I = 1; I <= min(N, M); I++) {
      A[I][I] = BETA;
    }
  } else {
    // Set the array to BETA on the diagonal and ALPHA on the
    // offdiagonal.

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= M; I++) {
        A[I][J] = ALPHA;
      }
    }
    for (I = 1; I <= min(M, N); I++) {
      A[I][I] = BETA;
    }
  }
}
