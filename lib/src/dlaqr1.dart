// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void dlaqr1(
  final int N,
  final Matrix<double> H_,
  final int LDH,
  final double SR1,
  final double SI1,
  final double SR2,
  final double SI2,
  final Array<double> V,
) {
  final H = H_.having(ld: LDH);
  const ZERO = 0.0;
  double H21S, H31S, S;

  // Quick return if possible

  if (N != 2 && N != 3) {
    return;
  }

  if (N == 2) {
    S = (H[1][1] - SR2).abs() + SI2.abs() + H[2][1].abs();
    if (S == ZERO) {
      V[1] = ZERO;
      V[2] = ZERO;
    } else {
      H21S = H[2][1] / S;
      V[1] = H21S * H[1][2] +
          (H[1][1] - SR1) * ((H[1][1] - SR2) / S) -
          SI1 * (SI2 / S);
      V[2] = H21S * (H[1][1] + H[2][2] - SR1 - SR2);
    }
  } else {
    S = (H[1][1] - SR2).abs() + SI2.abs() + H[2][1].abs() + H[3][1].abs();
    if (S == ZERO) {
      V[1] = ZERO;
      V[2] = ZERO;
      V[3] = ZERO;
    } else {
      H21S = H[2][1] / S;
      H31S = H[3][1] / S;
      V[1] = (H[1][1] - SR1) * ((H[1][1] - SR2) / S) -
          SI1 * (SI2 / S) +
          H[1][2] * H21S +
          H[1][3] * H31S;
      V[2] = H21S * (H[1][1] + H[2][2] - SR1 - SR2) + H[2][3] * H31S;
      V[3] = H31S * (H[1][1] + H[3][3] - SR1 - SR2) + H21S * H[3][2];
    }
  }
}
