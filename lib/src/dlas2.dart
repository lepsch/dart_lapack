// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';

void dlas2(
  final double F,
  final double G,
  final double H,
  final Box<double> SSMIN,
  final Box<double> SSMAX,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const ONE = 1.0;
  const TWO = 2.0;
  double AS, AT, AU, C, FA, FHMN, FHMX, GA, HA;

  FA = F.abs();
  GA = G.abs();
  HA = H.abs();
  FHMN = min(FA, HA);
  FHMX = max(FA, HA);
  if (FHMN == ZERO) {
    SSMIN.value = ZERO;
    if (FHMX == ZERO) {
      SSMAX.value = GA;
    } else {
      SSMAX.value =
          max(FHMX, GA) * sqrt(ONE + pow(min(FHMX, GA) / max(FHMX, GA), 2));
    }
  } else {
    if (GA < FHMX) {
      AS = ONE + FHMN / FHMX;
      AT = (FHMX - FHMN) / FHMX;
      AU = pow(GA / FHMX, 2).toDouble();
      C = TWO / (sqrt(AS * AS + AU) + sqrt(AT * AT + AU));
      SSMIN.value = FHMN * C;
      SSMAX.value = FHMX / C;
    } else {
      AU = FHMX / GA;
      if (AU == ZERO) {
        // Avoid possible harmful underflow if exponent range
        // asymmetric (true SSMIN may not underflow even if
        // AU underflows)

        SSMIN.value = (FHMN * FHMX) / GA;
        SSMAX.value = GA;
      } else {
        AS = ONE + FHMN / FHMX;
        AT = (FHMX - FHMN) / FHMX;
        C = ONE / (sqrt(ONE + pow(AS * AU, 2)) + sqrt(ONE + pow(AT * AU, 2)));
        SSMIN.value = (FHMN * C) * AU;
        SSMIN.value += SSMIN.value;
        SSMAX.value = GA / (C + C);
      }
    }
  }
}
