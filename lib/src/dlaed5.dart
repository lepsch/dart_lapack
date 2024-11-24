// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlaed5(
  final int I,
  final Array<double> D_,
  final Array<double> Z_,
  final Array<double> DELTA_,
  final double RHO,
  final Box<double> DLAM,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Z = Z_.having();
  final DELTA = DELTA_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, FOUR = 4.0;
  double B, C, DEL, TAU, TEMP, W;

  DEL = D[2] - D[1];
  if (I == 1) {
    W = ONE + TWO * RHO * (Z[2] * Z[2] - Z[1] * Z[1]) / DEL;
    if (W > ZERO) {
      B = DEL + RHO * (Z[1] * Z[1] + Z[2] * Z[2]);
      C = RHO * Z[1] * Z[1] * DEL;

      // B > ZERO, always

      TAU = TWO * C / (B + sqrt((B * B - FOUR * C).abs()));
      DLAM.value = D[1] + TAU;
      DELTA[1] = -Z[1] / TAU;
      DELTA[2] = Z[2] / (DEL - TAU);
    } else {
      B = -DEL + RHO * (Z[1] * Z[1] + Z[2] * Z[2]);
      C = RHO * Z[2] * Z[2] * DEL;
      if (B > ZERO) {
        TAU = -TWO * C / (B + sqrt(B * B + FOUR * C));
      } else {
        TAU = (B - sqrt(B * B + FOUR * C)) / TWO;
      }
      DLAM.value = D[2] + TAU;
      DELTA[1] = -Z[1] / (DEL + TAU);
      DELTA[2] = -Z[2] / TAU;
    }
    TEMP = sqrt(DELTA[1] * DELTA[1] + DELTA[2] * DELTA[2]);
    DELTA[1] /= TEMP;
    DELTA[2] /= TEMP;
  } else {
    // Now I=2

    B = -DEL + RHO * (Z[1] * Z[1] + Z[2] * Z[2]);
    C = RHO * Z[2] * Z[2] * DEL;
    if (B > ZERO) {
      TAU = (B + sqrt(B * B + FOUR * C)) / TWO;
    } else {
      TAU = TWO * C / (-B + sqrt(B * B + FOUR * C));
    }
    DLAM.value = D[2] + TAU;
    DELTA[1] = -Z[1] / (DEL + TAU);
    DELTA[2] = -Z[2] / TAU;
    TEMP = sqrt(DELTA[1] * DELTA[1] + DELTA[2] * DELTA[2]);
    DELTA[1] /= TEMP;
    DELTA[2] /= TEMP;
  }
}
