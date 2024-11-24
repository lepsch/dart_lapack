// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/nio.dart';

import '../test_driver.dart';

void chkxer(
  final String SRNAMT,
  final int INFOT,
  final Nout NOUT,
  final Box<bool> LERR,
  final Box<bool> OK, [
  final TestDriver? test,
]) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final reason =
      'Illegal value of parameter number ${INFOT.i2} not detected by ${SRNAMT.trim()}';
  test?.expect(LERR.value, true, reason: reason);
  if (!LERR.value) {
    NOUT.println(' *** $reason ***');
    OK.value = false;
  }
  LERR.value = false;
}
