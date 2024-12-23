// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlarz.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlatrz(
  final int M,
  final int N,
  final int L,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  int I;

  // Quick return if possible

  if (M == 0) {
    return;
  } else if (M == N) {
    for (I = 1; I <= N; I++) {
      TAU[I] = ZERO;
    }
    return;
  }

  for (I = M; I >= 1; I--) {
    // Generate elementary reflector H(i) to annihilate
    // [ A(i,i) A(i,n-l+1:n) ]

    dlarfg(L + 1, A.box(I, I), A(I, N - L + 1).asArray(), LDA, TAU.box(I));

    // Apply H(i) to A(1:i-1,i:n) from the right

    dlarz('Right', I - 1, N - I + 1, L, A(I, N - L + 1).asArray(), LDA, TAU[I],
        A(1, I), LDA, WORK);
  }
}
