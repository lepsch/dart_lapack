// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfg.dart';

void zgeql2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
    xerbla('ZGEQL2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = K; I >= 1; I--) {
    // Generate elementary reflector H(i) to annihilate
    // A(1:m-k+i-1,n-k+i)

    ALPHA.value = A[M - K + I][N - K + I];
    zlarfg(M - K + I, ALPHA, A(1, N - K + I).asArray(), 1, TAU(I));

    // Apply H(i)**H to A(1:m-k+i,1:n-k+i-1) from the left

    A[M - K + I][N - K + I] = Complex.one;
    zlarf('Left', M - K + I, N - K + I - 1, A(1, N - K + I).asArray(), 1,
        TAU[I].conjugate(), A, LDA, WORK);
    A[M - K + I][N - K + I] = ALPHA.value;
  }
}
