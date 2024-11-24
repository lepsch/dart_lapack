// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlasq6(
  final int I0,
  final int N0,
  final Array<double> Z_,
  final int PP,
  final Box<double> DMIN,
  final Box<double> DMIN1,
  final Box<double> DMIN2,
  final Box<double> DN,
  final Box<double> DNM1,
  final Box<double> DNM2,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having();
  const ZERO = 0.0;
  int J4, J4P2;
  double D, EMIN, SAFMIN, TEMP;

  if ((N0 - I0 - 1) <= 0) return;

  SAFMIN = dlamch('Safe minimum');
  J4 = 4 * I0 + PP - 3;
  EMIN = Z[J4 + 4];
  D = Z[J4];
  DMIN.value = D;

  if (PP == 0) {
    for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
      Z[J4 - 2] = D + Z[J4 - 1];
      if (Z[J4 - 2] == ZERO) {
        Z[J4] = ZERO;
        D = Z[J4 + 1];
        DMIN.value = D;
        EMIN = ZERO;
      } else if (SAFMIN * Z[J4 + 1] < Z[J4 - 2] &&
          SAFMIN * Z[J4 - 2] < Z[J4 + 1]) {
        TEMP = Z[J4 + 1] / Z[J4 - 2];
        Z[J4] = Z[J4 - 1] * TEMP;
        D *= TEMP;
      } else {
        Z[J4] = Z[J4 + 1] * (Z[J4 - 1] / Z[J4 - 2]);
        D = Z[J4 + 1] * (D / Z[J4 - 2]);
      }
      DMIN.value = min(DMIN.value, D);
      EMIN = min(EMIN, Z[J4]);
    }
  } else {
    for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
      Z[J4 - 3] = D + Z[J4];
      if (Z[J4 - 3] == ZERO) {
        Z[J4 - 1] = ZERO;
        D = Z[J4 + 2];
        DMIN.value = D;
        EMIN = ZERO;
      } else if (SAFMIN * Z[J4 + 2] < Z[J4 - 3] &&
          SAFMIN * Z[J4 - 3] < Z[J4 + 2]) {
        TEMP = Z[J4 + 2] / Z[J4 - 3];
        Z[J4 - 1] = Z[J4] * TEMP;
        D *= TEMP;
      } else {
        Z[J4 - 1] = Z[J4 + 2] * (Z[J4] / Z[J4 - 3]);
        D = Z[J4 + 2] * (D / Z[J4 - 3]);
      }
      DMIN.value = min(DMIN.value, D);
      EMIN = min(EMIN, Z[J4 - 1]);
    }
  }

  // Unroll last two steps.

  DNM2.value = D;
  DMIN2.value = DMIN.value;
  J4 = 4 * (N0 - 2) - PP;
  J4P2 = J4 + 2 * PP - 1;
  Z[J4 - 2] = DNM2.value + Z[J4P2];
  if (Z[J4 - 2] == ZERO) {
    Z[J4] = ZERO;
    DNM1.value = Z[J4P2 + 2];
    DMIN.value = DNM1.value;
    EMIN = ZERO;
  } else if (SAFMIN * Z[J4P2 + 2] < Z[J4 - 2] &&
      SAFMIN * Z[J4 - 2] < Z[J4P2 + 2]) {
    TEMP = Z[J4P2 + 2] / Z[J4 - 2];
    Z[J4] = Z[J4P2] * TEMP;
    DNM1.value = DNM2.value * TEMP;
  } else {
    Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
    DNM1.value = Z[J4P2 + 2] * (DNM2.value / Z[J4 - 2]);
  }
  DMIN.value = min(DMIN.value, DNM1.value);

  DMIN1.value = DMIN.value;
  J4 += 4;
  J4P2 = J4 + 2 * PP - 1;
  Z[J4 - 2] = DNM1.value + Z[J4P2];
  if (Z[J4 - 2] == ZERO) {
    Z[J4] = ZERO;
    DN.value = Z[J4P2 + 2];
    DMIN.value = DN.value;
    EMIN = ZERO;
  } else if (SAFMIN * Z[J4P2 + 2] < Z[J4 - 2] &&
      SAFMIN * Z[J4 - 2] < Z[J4P2 + 2]) {
    TEMP = Z[J4P2 + 2] / Z[J4 - 2];
    Z[J4] = Z[J4P2] * TEMP;
    DN.value = DNM1.value * TEMP;
  } else {
    Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
    DN.value = Z[J4P2 + 2] * (DNM1.value / Z[J4 - 2]);
  }
  DMIN.value = min(DMIN.value, DN.value);

  Z[J4 + 2] = DN.value;
  Z[4 * N0 - PP] = EMIN;
}
