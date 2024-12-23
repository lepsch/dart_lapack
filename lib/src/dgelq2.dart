// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgelq2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  int I, K;
  double AII;

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
    xerbla('DGELQ2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = 1; I <= K; I++) {
    // Generate elementary reflector H(i) to annihilate A(i,i+1:n)

    dlarfg(
        N - I + 1, A.box(I, I), A(I, min(I + 1, N)).asArray(), LDA, TAU.box(I));
    if (I < M) {
      // Apply H(i) to A(i+1:m,i:n) from the right

      AII = A[I][I];
      A[I][I] = ONE;
      dlarf('Right', M - I, N - I + 1, A(I, I).asArray(), LDA, TAU[I],
          A(I + 1, I), LDA, WORK);
      A[I][I] = AII;
    }
  }
}
