// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlarf.dart';

void zung2l(
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  final TAU = TAU_.having();
  int I, II, J, L;

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
    xerbla('ZUNG2L', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  // Initialise columns 1:n-k to columns of the unit matrix

  for (J = 1; J <= N - K; J++) {
    for (L = 1; L <= M; L++) {
      A[L][J] = Complex.zero;
    }
    A[M - N + J][J] = Complex.one;
  }

  for (I = 1; I <= K; I++) {
    II = N - K + I;

    // Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

    A[M - N + II][II] = Complex.one;
    zlarf('Left', M - N + II, II - 1, A(1, II).asArray(), 1, TAU[I], A, LDA,
        WORK);
    zscal(M - N + II - 1, -TAU[I], A(1, II).asArray(), 1);
    A[M - N + II][II] = Complex.one - TAU[I];

    // Set A(m-k+i+1:m,n-k+i) to zero

    for (L = M - N + II + 1; L <= M; L++) {
      A[L][II] = Complex.zero;
    }
  }
}
