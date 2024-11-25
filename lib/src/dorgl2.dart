// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorgl2(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, L;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < M) {
    INFO.value = -2;
  } else if (K < 0 || K > M) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DORGL2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M <= 0) return;

  if (K < M) {
    // Initialise rows k+1:m to rows of the unit matrix

    for (J = 1; J <= N; J++) {
      for (L = K + 1; L <= M; L++) {
        A[L][J] = ZERO;
      }
      if (J > K && J <= M) A[J][J] = ONE;
    }
  }

  for (I = K; I >= 1; I--) {
    // Apply H(i) to A(i:m,i:n) from the right

    if (I < N) {
      if (I < M) {
        A[I][I] = ONE;
        dlarf('Right', M - I, N - I + 1, A(I, I).asArray(), LDA, TAU[I],
            A(I + 1, I), LDA, WORK);
      }
      dscal(N - I, -TAU[I], A(I, I + 1).asArray(), LDA);
    }
    A[I][I] = ONE - TAU[I];

    // Set A(i,1:i-1) to zero

    for (L = 1; L <= I - 1; L++) {
      A[I][L] = ZERO;
    }
  }
}
