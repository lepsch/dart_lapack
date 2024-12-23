// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

double dsxt1(
  final int IJOB,
  final Array<double> D1_,
  final int N1,
  final Array<double> D2_,
  final int N2,
  final double ABSTOL,
  final double ULP,
  final double UNFL,
) {
  final D1 = D1_.having();
  final D2 = D2_.having();
  const ZERO = 0.0;
  int I, J;
  double TEMP1, TEMP2;

  TEMP1 = ZERO;

  J = 1;
  for (I = 1; I <= N1; I++) {
    while (D2[J] < D1[I] && J < N2) {
      J++;
    }
    if (J == 1) {
      TEMP2 = (D2[J] - D1[I]).abs();
      if (IJOB == 2) TEMP2 /= max(UNFL, ABSTOL + ULP * D1[I].abs());
    } else {
      TEMP2 = min((D2[J] - D1[I]).abs(), (D1[I] - D2[J - 1]).abs());
      if (IJOB == 2) TEMP2 /= max(UNFL, ABSTOL + ULP * D1[I].abs());
    }
    TEMP1 = max(TEMP1, TEMP2);
  }

  return TEMP1;
}
