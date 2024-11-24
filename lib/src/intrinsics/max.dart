// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math' as math;

T max<T extends num>(T a, T b) {
  if (a.isNaN) return b;
  if (b.isNaN) return a;
  return math.max(a, b);
}
