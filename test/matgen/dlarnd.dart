// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/matrix.dart';

import 'dlaran.dart';

double dlarnd(final int IDIST, final Array<int> ISEED_) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having();
  const ONE = 1.0, TWO = 2.0;
  const TWOPI = 6.28318530717958647692528676655900576839;
  double T1, T2;

  // Generate a real random number from a uniform (0,1) distribution

  T1 = dlaran(ISEED);

  if (IDIST == 1) {
    // uniform (0,1)

    return T1;
  } else if (IDIST == 2) {
    // uniform (-1,1)

    return TWO * T1 - ONE;
  } else if (IDIST == 3) {
    // normal (0,1)

    T2 = dlaran(ISEED);
    return sqrt(-TWO * log(T1)) * cos(TWOPI * T2);
  }
  return 0;
}
