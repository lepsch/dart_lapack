// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/dlamch.dart';

double dget06(final double RCOND, final double RCONDC) {
  const ZERO = 0.0, ONE = 1.0;

  final EPS = dlamch('Epsilon');
  if (RCOND > ZERO) {
    if (RCONDC > ZERO) {
      return max(RCOND, RCONDC) / min(RCOND, RCONDC) - (ONE - EPS);
    }

    return RCOND / EPS;
  }

  if (RCONDC > ZERO) {
    return RCONDC / EPS;
  }

  return ZERO;
}
