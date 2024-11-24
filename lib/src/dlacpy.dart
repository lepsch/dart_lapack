// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlacpy(
  final String UPLO,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= min(J, M); I++) {
        B[I][J] = A[I][J];
      }
    }
  } else if (lsame(UPLO, 'L')) {
    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= M; I++) {
        B[I][J] = A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= M; I++) {
        B[I][J] = A[I][J];
      }
    }
  }
}
