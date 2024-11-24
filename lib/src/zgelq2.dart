// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfg.dart';

void zgelq2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  int I, K;
  final ALPHA = Box(Complex.zero);

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZGELQ2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = 1; I <= K; I++) {
    // Generate elementary reflector H(i) to annihilate A(i,i+1:n)

    zlacgv(N - I + 1, A(I, I).asArray(), LDA);
    ALPHA.value = A[I][I];
    zlarfg(N - I + 1, ALPHA, A(I, min(I + 1, N)).asArray(), LDA, TAU(I));
    if (I < M) {
      // Apply H(i) to A(i+1:m,i:n) from the right

      A[I][I] = Complex.one;
      zlarf('Right', M - I, N - I + 1, A(I, I).asArray(), LDA, TAU[I],
          A(I + 1, I), LDA, WORK);
    }
    A[I][I] = ALPHA.value;
    zlacgv(N - I + 1, A(I, I).asArray(), LDA);
  }
}
