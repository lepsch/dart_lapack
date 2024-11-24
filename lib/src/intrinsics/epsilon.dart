// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/intrinsics/digits.dart';
import 'package:dart_lapack/src/intrinsics/radix.dart';

const double Function(double v) epsilon = epsilon64;

double epsilon32(final double v) => pow(radix32(v), 1 - digits32(v)).toDouble();
double epsilon64(final double v) => pow(radix64(v), 1 - digits64(v)).toDouble();
