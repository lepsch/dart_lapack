// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarra(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> E2_,
  final double SPLTOL,
  final double TNRM,
  final Box<int> NSPLIT,
  final Array<int> ISPLIT_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final E2 = E2_.having();
  final ISPLIT = ISPLIT_.having();
  const ZERO = 0.0;
  int I;
  double EABS, TMP1;

  INFO.value = 0;
  NSPLIT.value = 1;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  // Compute splitting points
  if (SPLTOL < ZERO) {
    // Criterion based on absolute off-diagonal value
    TMP1 = SPLTOL.abs() * TNRM;
    for (I = 1; I <= N - 1; I++) {
      EABS = E[I].abs();
      if (EABS <= TMP1) {
        E[I] = ZERO;
        E2[I] = ZERO;
        ISPLIT[NSPLIT.value] = I;
        NSPLIT.value++;
      }
    }
  } else {
    // Criterion that guarantees relative accuracy
    for (I = 1; I <= N - 1; I++) {
      EABS = E[I].abs();
      if (EABS <= SPLTOL * sqrt(D[I].abs()) * sqrt(D[I + 1].abs())) {
        E[I] = ZERO;
        E2[I] = ZERO;
        ISPLIT[NSPLIT.value] = I;
        NSPLIT.value++;
      }
    }
  }
  ISPLIT[NSPLIT.value] = N;
}
