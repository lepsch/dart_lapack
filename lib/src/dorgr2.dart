// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorgr2(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, II, J, L;

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
    xerbla('DORGR2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M <= 0) return;

  if (K < M) {
    // Initialise rows 1:m-k to rows of the unit matrix

    for (J = 1; J <= N; J++) {
      for (L = 1; L <= M - K; L++) {
        A[L][J] = ZERO;
      }
      if (J > N - M && J <= N - K) A[M - N + J][J] = ONE;
    }
  }

  for (I = 1; I <= K; I++) {
    II = M - K + I;

    // Apply H(i) to A(1:m-k+i,1:n-k+i) from the right

    A[II][N - M + II] = ONE;
    dlarf('Right', II - 1, N - M + II, A(II, 1).asArray(), LDA, TAU[I], A, LDA,
        WORK);
    dscal(N - M + II - 1, -TAU[I], A(II, 1).asArray(), LDA);
    A[II][N - M + II] = ONE - TAU[I];

    // Set A(m-k+i,n-k+i+1:n) to zero

    for (L = N - M + II + 1; L <= N; L++) {
      A[II][L] = ZERO;
    }
  }
}
