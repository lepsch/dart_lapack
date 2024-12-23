// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgeql2(
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
    xerbla('DGEQL2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = K; I >= 1; I--) {
    // Generate elementary reflector H(i) to annihilate
    // A(1:m-k+i-1,n-k+i)

    dlarfg(M - K + I, A.box(M - K + I, N - K + I), A(1, N - K + I).asArray(), 1,
        TAU.box(I));

    // Apply H(i) to A(1:m-k+i,1:n-k+i-1) from the left

    AII = A[M - K + I][N - K + I];
    A[M - K + I][N - K + I] = ONE;
    dlarf('Left', M - K + I, N - K + I - 1, A(1, N - K + I).asArray(), 1,
        TAU[I], A, LDA, WORK);
    A[M - K + I][N - K + I] = AII;
  }
}
