// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dger.dart';
import 'package:dart_lapack/src/blas/dtrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgeqrt2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  const ONE = 1.0e+00, ZERO = 0.0e+00;
  int I, K;
  double AII, ALPHA;

  // Test the input arguments

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -2;
  } else if (M < N) {
    INFO.value = -1;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DGEQRT2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = 1; I <= K; I++) {
    // Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)

    dlarfg(
        M - I + 1, A.box(I, I), A(min(I + 1, M), I).asArray(), 1, T.box(I, 1));
    if (I < N) {
      // Apply H(i) to A(I:M,I+1:N) from the left

      AII = A[I][I];
      A[I][I] = ONE;

      // W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]

      dgemv('T', M - I + 1, N - I, ONE, A(I, I + 1), LDA, A(I, I).asArray(), 1,
          ZERO, T(1, N).asArray(), 1);

      // A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H

      ALPHA = -T[I][1];
      dger(M - I + 1, N - I, ALPHA, A(I, I).asArray(), 1, T(1, N).asArray(), 1,
          A(I, I + 1), LDA);
      A[I][I] = AII;
    }
  }

  for (I = 2; I <= N; I++) {
    AII = A[I][I];
    A[I][I] = ONE;

    // T(1:I-1,I) := alpha * A(I:M,1:I-1)**T * A(I:M,I)

    ALPHA = -T[I][1];
    dgemv('T', M - I + 1, I - 1, ALPHA, A(I, 1), LDA, A(I, I).asArray(), 1,
        ZERO, T(1, I).asArray(), 1);
    A[I][I] = AII;

    // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

    dtrmv('U', 'N', 'N', I - 1, T, LDT, T(1, I).asArray(), 1);

    // T(I,I) = tau(I)

    T[I][I] = T[I][1];
    T[I][1] = ZERO;
  }
}
