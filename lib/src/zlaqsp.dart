// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlaqsp(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<double> S_,
  final double SCOND,
  final double AMAX,
  final Box<String> EQUED,
) {
  final AP = AP_.having();
  final S = S_.having();
  const ONE = 1.0, THRESH = 0.1;
  int I, J, JC;
  double CJ, LARGE, SMALL;

  // Quick return if possible

  if (N <= 0) {
    EQUED.value = 'N';
    return;
  }

  // Initialize LARGE and SMALL.

  SMALL = dlamch('Safe minimum') / dlamch('Precision');
  LARGE = ONE / SMALL;

  if (SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE) {
    // No equilibration

    EQUED.value = 'N';
  } else {
    // Replace A by diag(S) * A * diag(S).

    if (lsame(UPLO, 'U')) {
      // Upper triangle of A is stored.

      JC = 1;
      for (J = 1; J <= N; J++) {
        CJ = S[J];
        for (I = 1; I <= J; I++) {
          AP[JC + I - 1] = (CJ * S[I]).toComplex() * AP[JC + I - 1];
        }
        JC += J;
      }
    } else {
      // Lower triangle of A is stored.

      JC = 1;
      for (J = 1; J <= N; J++) {
        CJ = S[J];
        for (I = J; I <= N; I++) {
          AP[JC + I - J] = (CJ * S[I]).toComplex() * AP[JC + I - J];
        }
        JC += N - J + 1;
      }
    }
    EQUED.value = 'Y';
  }
}
