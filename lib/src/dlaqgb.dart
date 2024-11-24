// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaqgb(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> R_,
  final Array<double> C_,
  final double ROWCND,
  final double COLCND,
  final double AMAX,
  final Box<String> EQUED,
) {
  final AB = AB_.having(ld: LDAB);
  final C = C_.having();
  final R = R_.having();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, THRESH = 0.1;
  int I, J;
  double CJ, LARGE, SMALL;

  // Quick return if possible

  if (M <= 0 || N <= 0) {
    EQUED.value = 'N';
    return;
  }

  // Initialize LARGE and SMALL.

  SMALL = dlamch('Safe minimum') / dlamch('Precision');
  LARGE = ONE / SMALL;

  if (ROWCND >= THRESH && AMAX >= SMALL && AMAX <= LARGE) {
    // No row scaling

    if (COLCND >= THRESH) {
      // No column scaling

      EQUED.value = 'N';
    } else {
      // Column scaling

      for (J = 1; J <= N; J++) {
        CJ = C[J];
        for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
          AB[KU + 1 + I - J][J] = CJ * AB[KU + 1 + I - J][J];
        }
      }
      EQUED.value = 'C';
    }
  } else if (COLCND >= THRESH) {
    // Row scaling, no column scaling

    for (J = 1; J <= N; J++) {
      for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
        AB[KU + 1 + I - J][J] = R[I] * AB[KU + 1 + I - J][J];
      }
    }
    EQUED.value = 'R';
  } else {
    // Row and column scaling

    for (J = 1; J <= N; J++) {
      CJ = C[J];
      for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
        AB[KU + 1 + I - J][J] = CJ * R[I] * AB[KU + 1 + I - J][J];
      }
    }
    EQUED.value = 'B';
  }
}
