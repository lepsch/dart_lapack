// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgttrf(
  final int N,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DU2_,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final IPIV = IPIV_.having();
  const ZERO = 0.0;
  int I;
  Complex FACT, TEMP;

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('ZGTTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Initialize IPIV(i) = i and DU2(i) = 0

  for (I = 1; I <= N; I++) {
    IPIV[I] = I;
  }
  for (I = 1; I <= N - 2; I++) {
    DU2[I] = Complex.zero;
  }

  for (I = 1; I <= N - 2; I++) {
    if (D[I].cabs1() >= DL[I].cabs1()) {
      // No row interchange required, eliminate DL(I)

      if (D[I].cabs1() != ZERO) {
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
    if (D[I].cabs1() >= DL[I].cabs1()) {
      if (D[I].cabs1() != ZERO) {
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
    if (D[I].cabs1() == ZERO) {
      INFO.value = I;
      return;
    }
  }
}
