// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/intrinsics.dart';
import 'package:dart_lapack/src/box.dart';

void dlabad(Box<double> SMALL, Box<double> LARGE) {
  // If it looks like we're on a Cray, take the square root of
  // SMALL and LARGE to avoid overflow and underflow problems.

  if (log10(LARGE.value) > 2000.0) {
    SMALL.value = sqrt(SMALL.value);
    LARGE.value = sqrt(LARGE.value);
  }
}
