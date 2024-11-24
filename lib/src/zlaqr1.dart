// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlaqr1(
  final int N,
  final Matrix<Complex> H_,
  final int LDH,
  final Complex S1,
  final Complex S2,
  final Array<Complex> V_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final V = V_.having();
  const RZERO = 0.0;
  Complex H21S, H31S;
  double S;

  // Quick return if possible

  if (N != 2 && N != 3) {
    return;
  }

  if (N == 2) {
    S = (H[1][1] - S2).cabs1() + H[2][1].cabs1();
    if (S == RZERO) {
      V[1] = Complex.zero;
      V[2] = Complex.zero;
    } else {
      H21S = H[2][1] / S.toComplex();
      V[1] = H21S * H[1][2] + (H[1][1] - S1) * ((H[1][1] - S2) / S.toComplex());
      V[2] = H21S * (H[1][1] + H[2][2] - S1 - S2);
    }
  } else {
    S = (H[1][1] - S2).cabs1() + H[2][1].cabs1() + H[3][1].cabs1();
    if (S == RZERO) {
      V[1] = Complex.zero;
      V[2] = Complex.zero;
      V[3] = Complex.zero;
    } else {
      H21S = H[2][1] / S.toComplex();
      H31S = H[3][1] / S.toComplex();
      V[1] = (H[1][1] - S1) * ((H[1][1] - S2) / S.toComplex()) +
          H[1][2] * H21S +
          H[1][3] * H31S;
      V[2] = H21S * (H[1][1] + H[2][2] - S1 - S2) + H[2][3] * H31S;
      V[3] = H31S * (H[1][1] + H[3][3] - S1 - S2) + H21S * H[3][2];
    }
  }
}
