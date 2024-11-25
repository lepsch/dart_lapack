// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dpoequb(
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> S_,
  final Box<double> SCOND,
  final Box<double> AMAX,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final S = S_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I;
  double SMIN, BASE, TMP;

  // Test the input parameters.

  // Positive definite only performs 1 pass of equilibration.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (LDA < max(1, N)) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('DPOEQUB', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) {
    SCOND.value = ONE;
    AMAX.value = ZERO;
    return;
  }

  BASE = dlamch('B');
  TMP = -0.5 / log(BASE);

  // Find the minimum and maximum diagonal elements.

  S[1] = A[1][1];
  SMIN = S[1];
  AMAX.value = S[1];
  for (I = 2; I <= N; I++) {
    S[I] = A[I][I];
    SMIN = min(SMIN, S[I]);
    AMAX.value = max(AMAX.value, S[I]);
  }

  if (SMIN <= ZERO) {
    // Find the first non-positive diagonal element and return.

    for (I = 1; I <= N; I++) {
      if (S[I] <= ZERO) {
        INFO.value = I;
        return;
      }
    }
  } else {
    // Set the scale factors to the reciprocals
    // of the diagonal elements.

    for (I = 1; I <= N; I++) {
      S[I] = pow(BASE, (TMP * log(S[I])).toInt()).toDouble();
    }

    // Compute SCOND = min(S(I)) / max(S(I)).

    SCOND.value = sqrt(SMIN) / sqrt(AMAX.value);
  }
}
