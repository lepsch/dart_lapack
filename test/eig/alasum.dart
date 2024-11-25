// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/nio.dart';

void alasum(
  final String TYPE,
  final Nout NOUT,
  final int NFAIL,
  final int NRUN,
  final int NERRS,
) {
  if (NFAIL > 0) {
    NOUT.println(
        ' ${TYPE.a3}: ${NFAIL.i6} out of ${NRUN.i6} tests failed to pass the threshold');
  } else {
    NOUT.println(
        '\n All tests for ${TYPE.a3} routines passed the threshold ( ${NRUN.i6} tests run)');
  }
  if (NERRS > 0) {
    NOUT.println('${' ' * 6}${NERRS.i6} error messages recorded');
  }
}
