// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/intrinsics/sign.dart';

bool dlctes(final double ZR, final double ZI, final double D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;

  if (D == ZERO) {
    return (ZR < ZERO);
  } else {
    return (sign(ONE, ZR) != sign(ONE, D));
  }
}
