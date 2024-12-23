// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math' as math;

num pow(num x, num exponent) {
  final e = exponent.truncateToDouble();
  if (e == exponent) {
    return _powi(x, e.toInt());
  }
  return math.pow(x, exponent);
}

num _powi(num a, int b) {
  num pow, x;
  int n;
  int u;

  n = b;
  x = a;
  pow = 1;
  if (n != 0) {
    if (n < 0) {
      u = -n;
      x = pow / x;
    } else {
      u = n;
    }
    for (;;) {
      if (u & 1 != 0) pow *= x;
      u >>>= 1;
      if (u != 0) {
        x *= x;
      } else {
        break;
      }
    }
  }
  return pow;
}
