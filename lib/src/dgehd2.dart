// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgehd2(
  final int N,
  final int ILO,
  final int IHI,
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
  int I;
  double AII;

  // Test the input parameters

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -2;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DGEHD2', -INFO.value);
    return;
  }

  for (I = ILO; I <= IHI - 1; I++) {
    // Compute elementary reflector H(i) to annihilate A[i+2:ihi][i]

    dlarfg(
        IHI - I, A.box(I + 1, I), A(min(I + 2, N), I).asArray(), 1, TAU.box(I));
    AII = A[I + 1][I];
    A[I + 1][I] = ONE;

    // Apply H(i) to A[1:ihi][i+1:ihi] from the right

    dlarf('Right', IHI, IHI - I, A(I + 1, I).asArray(), 1, TAU[I], A(1, I + 1),
        LDA, WORK);

    // Apply H(i) to A[i+1:ihi][i+1:n] from the left

    dlarf('Left', IHI - I, N - I, A(I + 1, I).asArray(), 1, TAU[I],
        A(I + 1, I + 1), LDA, WORK);

    A[I + 1][I] = AII;
  }
}
