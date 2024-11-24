// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/intrinsics/epsilon.dart';

double sroundup_lwork(final int LWORK) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  var result = LWORK.toDouble();

  if (result.toInt() < LWORK) {
    // Force round up of LWORK
    result *= (1.0 + epsilon(0.0));
  }
  return result;
}
