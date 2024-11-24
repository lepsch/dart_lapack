// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/dlapy2.dart';

import 'common.dart';

bool dslect(final double ZR, final double ZI) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I;
  double RMIN, X;
  const ZERO = 0.0;

  if (sslct.SELOPT == 0) {
    return (ZR < ZERO);
  }
  RMIN = dlapy2(ZR - sslct.SELWR[1], ZI - sslct.SELWI[1]);
  var result = sslct.SELVAL[1];
  for (I = 2; I <= sslct.SELDIM; I++) {
    X = dlapy2(ZR - sslct.SELWR[I], ZI - sslct.SELWI[I]);
    if (X <= RMIN) {
      RMIN = X;
      result = sslct.SELVAL[I];
    }
  }

  return result;
}
