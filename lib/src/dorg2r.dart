// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorg2r(
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
  int I, J, L;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N > M) {
    INFO.value = -2;
  } else if (K < 0 || K > N) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DORG2R', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  // Initialise columns k+1:n to columns of the unit matrix

  for (J = K + 1; J <= N; J++) {
    for (L = 1; L <= M; L++) {
      A[L][J] = ZERO;
    }
    A[J][J] = ONE;
  }

  for (I = K; I >= 1; I--) {
    // Apply H(i) to A(i:m,i:n) from the left

    if (I < N) {
      A[I][I] = ONE;
      dlarf('Left', M - I + 1, N - I, A(I, I).asArray(), 1, TAU[I], A(I, I + 1),
          LDA, WORK);
    }
    if (I < M) dscal(M - I, -TAU[I], A(I + 1, I).asArray(), 1);
    A[I][I] = ONE - TAU[I];

    // Set A(1:i-1,i) to zero

    for (L = 1; L <= I - 1; L++) {
      A[L][I] = ZERO;
    }
  }
}
