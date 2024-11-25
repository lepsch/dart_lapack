// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgttrf(
  final int N,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final IPIV = IPIV_.having();
  const ZERO = 0.0;
  int I;
  double FACT, TEMP;

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('DGTTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Initialize IPIV(i) = i and DU2(I) = 0

  for (I = 1; I <= N; I++) {
    IPIV[I] = I;
  }
  for (I = 1; I <= N - 2; I++) {
    DU2[I] = ZERO;
  }

  for (I = 1; I <= N - 2; I++) {
    if (D[I].abs() >= DL[I].abs()) {
      // No row interchange required, eliminate DL(I)

      if (D[I] != ZERO) {
        FACT = DL[I] / D[I];
        DL[I] = FACT;
        D[I + 1] -= FACT * DU[I];
      }
    } else {
      // Interchange rows I and I+1, eliminate DL(I)

      FACT = D[I] / DL[I];
      D[I] = DL[I];
      DL[I] = FACT;
      TEMP = DU[I];
      DU[I] = D[I + 1];
      D[I + 1] = TEMP - FACT * D[I + 1];
      DU2[I] = DU[I + 1];
      DU[I + 1] = -FACT * DU[I + 1];
      IPIV[I] = I + 1;
    }
  }
  if (N > 1) {
    I = N - 1;
    if (D[I].abs() >= DL[I].abs()) {
      if (D[I] != ZERO) {
        FACT = DL[I] / D[I];
        DL[I] = FACT;
        D[I + 1] -= FACT * DU[I];
      }
    } else {
      FACT = D[I] / DL[I];
      D[I] = DL[I];
      DL[I] = FACT;
      TEMP = DU[I];
      DU[I] = D[I + 1];
      D[I + 1] = TEMP - FACT * D[I + 1];
      IPIV[I] = I + 1;
    }
  }

  // Check for a zero on the diagonal of U.

  for (I = 1; I <= N; I++) {
    if (D[I] == ZERO) {
      INFO.value = I;
      break;
    }
  }
}
