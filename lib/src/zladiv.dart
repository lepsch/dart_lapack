// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dladiv.dart';

Complex zladiv(final Complex X, final Complex Y) {
  final ZI = Box(0.0), ZR = Box(0.0);
  dladiv(X.real, X.imaginary, Y.real, Y.imaginary, ZR, ZI);
  return Complex(ZR.value, ZI.value);
}
