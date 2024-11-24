// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'dlaran.dart';

Complex zlarnd(final int IDIST, final Array<int> ISEED) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const TWOPI = 6.28318530717958647692528676655900576839;
  double T1, T2;

  // Generate a pair of real random numbers from a uniform (0,1)
  // distribution

  T1 = dlaran(ISEED);
  T2 = dlaran(ISEED);

  if (IDIST == 1) {
    // real and imaginary parts each uniform (0,1)

    return Complex(T1, T2);
  } else if (IDIST == 2) {
    // real and imaginary parts each uniform (-1,1)

    return Complex(TWO * T1 - ONE, TWO * T2 - ONE);
  } else if (IDIST == 3) {
    // real and imaginary parts each normal (0,1)

    return sqrt(-TWO * log(T1)).toComplex() * Complex(ZERO, TWOPI * T2).exp();
  } else if (IDIST == 4) {
    // uniform distribution on the unit disc abs(z) <= 1

    return sqrt(T1).toComplex() * Complex(ZERO, TWOPI * T2).exp();
  } else if (IDIST == 5) {
    // uniform distribution on the unit circle abs(z) = 1

    return Complex(ZERO, TWOPI * T2).exp();
  }
  assert(false);
  return Complex.zero;
}
