// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/dlamch.dart';

double dlarmm(
  final double ANORM,
  final double BNORM,
  final double CNORM,
) {
  const ONE = 1.0, HALF = 0.5, FOUR = 4.0;
  double BIGNUM, SMLNUM;

  // Determine machine dependent parameters to control overflow.

  SMLNUM = dlamch('Safe minimum') / dlamch('Precision');
  BIGNUM = (ONE / SMLNUM) / FOUR;

  // Compute a scale factor.

  if (BNORM <= ONE) {
    if (ANORM * BNORM > BIGNUM - CNORM) {
      return HALF;
    }
  } else {
    if (ANORM > (BIGNUM - CNORM) / BNORM) {
      return HALF / BNORM;
    }
  }
  return ONE;
}
