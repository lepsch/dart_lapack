// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/matrix.dart';

double dla_gbrpvgrw(
  final int N,
  final int KL,
  final int KU,
  final int NCOLS,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> AFB_,
  final int LDAFB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
  int I, J, KD;
  double AMAX, UMAX, RPVGRW;

  RPVGRW = 1.0;

  KD = KU + 1;
  for (J = 1; J <= NCOLS; J++) {
    AMAX = 0.0;
    UMAX = 0.0;
    for (I = max(J - KU, 1); I <= min(J + KL, N); I++) {
      AMAX = max(AB[KD + I - J][J].abs(), AMAX);
    }
    for (I = max(J - KU, 1); I <= J; I++) {
      UMAX = max(AFB[KD + I - J][J].abs(), UMAX);
    }
    if (UMAX != 0.0) {
      RPVGRW = min(AMAX / UMAX, RPVGRW);
    }
  }
  return RPVGRW;
}
