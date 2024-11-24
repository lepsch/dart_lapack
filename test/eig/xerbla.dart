// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/format_specifiers_extensions.dart';

import '../test_driver.dart';
import 'common.dart';

void Function(String SRNAME, int INFO) xerbla(final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  return (final String SRNAME, final int INFO) {
    infoc.LERR.value = true;

    var reason = infoc.INFOT != 0
        ? 'XERBLA was called from ${srnamc.SRNAMT} with INFO = ${INFO.i6} instead of ${infoc.INFOT.i2}'
        : 'On entry to $SRNAME parameter number ${INFO.i6} had an illegal value';
    test.expect(INFO, infoc.INFOT, reason: reason);
    if (INFO != infoc.INFOT) {
      infoc.NOUT.println(' *** $reason ***');
      infoc.OK.value = false;
    }
    reason =
        'XERBLA was called with SRNAME = $SRNAME instead of ${srnamc.SRNAMT.a6}';
    test.expect(SRNAME, srnamc.SRNAMT, reason: reason);
    if (SRNAME != srnamc.SRNAMT) {
      infoc.NOUT.println(' *** $reason ***');
      infoc.OK.value = false;
    }
  };
}
