// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/complex.dart';

import 'common.dart';

bool zslect(final Complex Z) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const ZERO = 0.0;

  if (sslct.SELOPT == 0) {
    return Z.real < ZERO;
  }

  var RMIN = (Z - Complex(sslct.SELWR[1], sslct.SELWI[1])).abs();
  var result = sslct.SELVAL[1];
  for (var I = 2; I <= sslct.SELDIM; I++) {
    final X = (Z - Complex(sslct.SELWR[I], sslct.SELWI[I])).abs();
    if (X <= RMIN) {
      RMIN = X;
      result = sslct.SELVAL[I];
    }
  }
  return result;
}
