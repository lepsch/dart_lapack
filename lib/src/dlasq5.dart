// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void dlasq5(
  final int I0,
  final int N0,
  final Array<double> Z_,
  final int PP,
  final Box<double> TAU,
  final double SIGMA,
  final Box<double> DMIN,
  final Box<double> DMIN1,
  final Box<double> DMIN2,
  final Box<double> DN,
  final Box<double> DNM1,
  final Box<double> DNM2,
  final bool IEEE,
  final double EPS,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having();
  const ZERO = 0.0, HALF = 0.5;
  int J4, J4P2;
  double D, EMIN, TEMP, DTHRESH;

  if ((N0 - I0 - 1) <= 0) return;

  DTHRESH = EPS * (SIGMA + TAU.value);
  if (TAU.value < DTHRESH * HALF) TAU.value = ZERO;
  if (TAU.value != ZERO) {
    J4 = 4 * I0 + PP - 3;
    EMIN = Z[J4 + 4];
    D = Z[J4] - TAU.value;
    DMIN.value = D;
    DMIN1.value = -Z[J4];

    if (IEEE) {
      // Code for IEEE arithmetic.

      if (PP == 0) {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 2] = D + Z[J4 - 1];
          TEMP = Z[J4 + 1] / Z[J4 - 2];
          D = D * TEMP - TAU.value;
          DMIN.value = min(DMIN.value, D);
          Z[J4] = Z[J4 - 1] * TEMP;
          EMIN = min(Z[J4], EMIN);
        }
      } else {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 3] = D + Z[J4];
          TEMP = Z[J4 + 2] / Z[J4 - 3];
          D = D * TEMP - TAU.value;
          DMIN.value = min(DMIN.value, D);
          Z[J4 - 1] = Z[J4] * TEMP;
          EMIN = min(Z[J4 - 1], EMIN);
        }
      }

      // Unroll last two steps.

      DNM2.value = D;
      DMIN2.value = DMIN.value;
      J4 = 4 * (N0 - 2) - PP;
      J4P2 = J4 + 2 * PP - 1;
      Z[J4 - 2] = DNM2.value + Z[J4P2];
      Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
      DNM1.value = Z[J4P2 + 2] * (DNM2.value / Z[J4 - 2]) - TAU.value;
      DMIN.value = min(DMIN.value, DNM1.value);

      DMIN1.value = DMIN.value;
      J4 += 4;
      J4P2 = J4 + 2 * PP - 1;
      Z[J4 - 2] = DNM1.value + Z[J4P2];
      Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
      DN.value = Z[J4P2 + 2] * (DNM1.value / Z[J4 - 2]) - TAU.value;
      DMIN.value = min(DMIN.value, DN.value);
    } else {
      // Code for non IEEE arithmetic.

      if (PP == 0) {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 2] = D + Z[J4 - 1];
          if (D < ZERO) {
            return;
          } else {
            Z[J4] = Z[J4 + 1] * (Z[J4 - 1] / Z[J4 - 2]);
            D = Z[J4 + 1] * (D / Z[J4 - 2]) - TAU.value;
          }
          DMIN.value = min(DMIN.value, D);
          EMIN = min(EMIN, Z[J4]);
        }
      } else {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 3] = D + Z[J4];
          if (D < ZERO) {
            return;
          } else {
            Z[J4 - 1] = Z[J4 + 2] * (Z[J4] / Z[J4 - 3]);
            D = Z[J4 + 2] * (D / Z[J4 - 3]) - TAU.value;
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
      if (DNM2.value < ZERO) {
        return;
      } else {
        Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
        DNM1.value = Z[J4P2 + 2] * (DNM2.value / Z[J4 - 2]) - TAU.value;
      }
      DMIN.value = min(DMIN.value, DNM1.value);

      DMIN1.value = DMIN.value;
      J4 += 4;
      J4P2 = J4 + 2 * PP - 1;
      Z[J4 - 2] = DNM1.value + Z[J4P2];
      if (DNM1.value < ZERO) {
        return;
      } else {
        Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
        DN.value = Z[J4P2 + 2] * (DNM1.value / Z[J4 - 2]) - TAU.value;
      }
      DMIN.value = min(DMIN.value, DN.value);
    }
  } else {
    // This is the version that sets d's to zero if they are small enough
    J4 = 4 * I0 + PP - 3;
    EMIN = Z[J4 + 4];
    D = Z[J4] - TAU.value;
    DMIN.value = D;
    DMIN1.value = -Z[J4];
    if (IEEE) {
      // Code for IEEE arithmetic.

      if (PP == 0) {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 2] = D + Z[J4 - 1];
          TEMP = Z[J4 + 1] / Z[J4 - 2];
          D = D * TEMP - TAU.value;
          if (D < DTHRESH) D = ZERO;
          DMIN.value = min(DMIN.value, D);
          Z[J4] = Z[J4 - 1] * TEMP;
          EMIN = min(Z[J4], EMIN);
        }
      } else {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 3] = D + Z[J4];
          TEMP = Z[J4 + 2] / Z[J4 - 3];
          D = D * TEMP - TAU.value;
          if (D < DTHRESH) D = ZERO;
          DMIN.value = min(DMIN.value, D);
          Z[J4 - 1] = Z[J4] * TEMP;
          EMIN = min(Z[J4 - 1], EMIN);
        }
      }

      // Unroll last two steps.

      DNM2.value = D;
      DMIN2.value = DMIN.value;
      J4 = 4 * (N0 - 2) - PP;
      J4P2 = J4 + 2 * PP - 1;
      Z[J4 - 2] = DNM2.value + Z[J4P2];
      Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
      DNM1.value = Z[J4P2 + 2] * (DNM2.value / Z[J4 - 2]) - TAU.value;
      DMIN.value = min(DMIN.value, DNM1.value);

      DMIN1.value = DMIN.value;
      J4 += 4;
      J4P2 = J4 + 2 * PP - 1;
      Z[J4 - 2] = DNM1.value + Z[J4P2];
      Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
      DN.value = Z[J4P2 + 2] * (DNM1.value / Z[J4 - 2]) - TAU.value;
      DMIN.value = min(DMIN.value, DN.value);
    } else {
      // Code for non IEEE arithmetic.

      if (PP == 0) {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 2] = D + Z[J4 - 1];
          if (D < ZERO) {
            return;
          } else {
            Z[J4] = Z[J4 + 1] * (Z[J4 - 1] / Z[J4 - 2]);
            D = Z[J4 + 1] * (D / Z[J4 - 2]) - TAU.value;
          }
          if (D < DTHRESH) D = ZERO;
          DMIN.value = min(DMIN.value, D);
          EMIN = min(EMIN, Z[J4]);
        }
      } else {
        for (J4 = 4 * I0; J4 <= 4 * (N0 - 3); J4 += 4) {
          Z[J4 - 3] = D + Z[J4];
          if (D < ZERO) {
            return;
          } else {
            Z[J4 - 1] = Z[J4 + 2] * (Z[J4] / Z[J4 - 3]);
            D = Z[J4 + 2] * (D / Z[J4 - 3]) - TAU.value;
          }
          if (D < DTHRESH) D = ZERO;
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
      if (DNM2.value < ZERO) {
        return;
      } else {
        Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
        DNM1.value = Z[J4P2 + 2] * (DNM2.value / Z[J4 - 2]) - TAU.value;
      }
      DMIN.value = min(DMIN.value, DNM1.value);

      DMIN1.value = DMIN.value;
      J4 += 4;
      J4P2 = J4 + 2 * PP - 1;
      Z[J4 - 2] = DNM1.value + Z[J4P2];
      if (DNM1.value < ZERO) {
        return;
      } else {
        Z[J4] = Z[J4P2 + 2] * (Z[J4P2] / Z[J4 - 2]);
        DN.value = Z[J4P2 + 2] * (DNM1.value / Z[J4 - 2]) - TAU.value;
      }
      DMIN.value = min(DMIN.value, DN.value);
    }
  }

  Z[J4 + 2] = DN.value;
  Z[4 * N0 - PP] = EMIN;
}
