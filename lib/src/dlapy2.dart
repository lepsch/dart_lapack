// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/install/dlamch.dart';

double dlapy2(final double X, final double Y) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const zero = 0.0;
  const ONE = 1.0;

  var result = 0.0;
  final X_IS_NAN = disnan(X);
  final Y_IS_NAN = disnan(Y);
  if (X_IS_NAN) result = X;
  if (Y_IS_NAN) result = Y;
  final HUGEVAL = dlamch('Overflow');

  if (!(X_IS_NAN || Y_IS_NAN)) {
    final xabs = X.abs();
    final yabs = Y.abs();
    final W = max(xabs, yabs);
    final Z = min(xabs, yabs);
    if (Z == zero || W > HUGEVAL) {
      result = W;
    } else {
      result = W * sqrt(ONE + pow(Z / W, 2));
    }
  }
  return result;
}
