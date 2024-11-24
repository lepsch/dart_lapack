// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:typed_data';

import 'package:dart_lapack/src/intrinsics/digits.dart';
import 'package:dart_lapack/src/intrinsics/maxexponent.dart';

const int Function(double x) exponent = exponent64;

int exponent64(final double x) {
  final bias = maxexponent(x) - 1;
  final d = digits(x) - 1;
  return ((Float64List.fromList([x]).buffer.asUint32List(0).last &
              0x7fffffff) >>>
          (d - 32)) -
      bias +
      1;
}

int exponent32(final double x) {
  final bias = maxexponent32(x) - 1;
  final d = digits32(x) - 1;
  return ((Float64List.fromList([x]).buffer.asUint32List(0).last &
              0x7fffffff) >>>
          (d - 32)) -
      bias +
      1;
}
