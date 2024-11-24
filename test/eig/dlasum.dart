// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/nio.dart';

void dlasum(
  final String TYPE,
  final Nout IOUNIT,
  final int IE,
  final int NRUN,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  if (IE > 0) {
    IOUNIT.println(
        ' ${TYPE.a3}: ${IE.i4} out of ${NRUN.i5} tests failed to pass the threshold');
  } else {
    IOUNIT.println(
        '\n All tests for ${TYPE.a3} passed the threshold ( ${NRUN.i5} tests run)');
  }
}
