// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/ddot.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlauu2(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  const ONE = 1.0;
  bool UPPER;
  int I;
  double AII;

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
    xerbla('DLAUU2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Compute the product U * U**T.

    for (I = 1; I <= N; I++) {
      AII = A[I][I];
      if (I < N) {
        A[I][I] =
            ddot(N - I + 1, A(I, I).asArray(), LDA, A(I, I).asArray(), LDA);
        dgemv('No transpose', I - 1, N - I, ONE, A(1, I + 1), LDA,
            A(I, I + 1).asArray(), LDA, AII, A(1, I).asArray(), 1);
      } else {
        dscal(I, AII, A(1, I).asArray(), 1);
      }
    }
  } else {
    // Compute the product L**T * L.

    for (I = 1; I <= N; I++) {
      AII = A[I][I];
      if (I < N) {
        A[I][I] = ddot(N - I + 1, A(I, I).asArray(), 1, A(I, I).asArray(), 1);
        dgemv('Transpose', N - I, I - 1, ONE, A(I + 1, 1), LDA,
            A(I + 1, I).asArray(), 1, AII, A(I, 1).asArray(), LDA);
      } else {
        dscal(I, AII, A(I, 1).asArray(), LDA);
      }
    }
  }
}
