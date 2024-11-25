// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';

bool zlctes(final Complex Z, final Complex D) {
  const ZERO = 0.0, ONE = 1.0;

  if (D == Complex.zero) {
    return Z.real < ZERO;
  }

  if (Z.real == ZERO || D.real == ZERO) {
    return sign(ONE, Z.imaginary) != sign(ONE, D.imaginary);
  }

  if (Z.imaginary == ZERO || D.imaginary == ZERO) {
    return sign(ONE, Z.real) != sign(ONE, D.real);
  }

  final ZMAX = max((Z.real).abs(), Z.imaginary.abs());
  return (Z.real / ZMAX) * D.real + (Z.imaginary / ZMAX) * D.imaginary < ZERO;
}
