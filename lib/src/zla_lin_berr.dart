// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zla_lin_berr(
  final int N,
  final int NZ,
  final int NRHS,
  final Matrix<Complex> RES_,
  final Matrix<double> AYB_,
  final Array<double> BERR_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RES = RES_.having(ld: N);
  final AYB = AYB_.having(ld: N);
  final BERR = BERR_.having(length: NRHS);
  double TMP, SAFE1;
  int I, J;

  // Adding SAFE1 to the numerator guards against spuriously zero
  // residuals.  A similar safeguard is in the CLA_yyAMV routine used
  // to compute AYB.

  SAFE1 = dlamch('Safe minimum');
  SAFE1 = (NZ + 1) * SAFE1;

  for (J = 1; J <= NRHS; J++) {
    BERR[J] = 0.0;
    for (I = 1; I <= N; I++) {
      if (AYB[I][J] != 0.0) {
        TMP = (SAFE1 + RES[I][J].cabs1()) / AYB[I][J];
        BERR[J] = max(BERR[J], TMP);
      }

      // If AYB is exactly 0.0 (and if computed by CLA_yyAMV), then we know
      // the true residual also must be exactly 0.0.
    }
  }
}
