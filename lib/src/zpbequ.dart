// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zpbequ(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> S_,
  final Box<double> SCOND,
  final Box<double> AMAX,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final S = S_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool UPPER;
  int I, J;
  double SMIN;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (LDAB < KD + 1) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZPBEQU', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) {
    SCOND.value = ONE;
    AMAX.value = ZERO;
    return;
  }

  if (UPPER) {
    J = KD + 1;
  } else {
    J = 1;
  }

  // Initialize SMIN and AMAX.

  S[1] = AB[J][1].real;
  SMIN = S[1];
  AMAX.value = S[1];

  // Find the minimum and maximum diagonal elements.

  for (I = 2; I <= N; I++) {
    S[I] = AB[J][I].real;
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
      S[I] = ONE / sqrt(S[I]);
    }

    // Compute SCOND = min(S(I)) / max(S(I))

    SCOND.value = sqrt(SMIN) / sqrt(AMAX.value);
  }
}
