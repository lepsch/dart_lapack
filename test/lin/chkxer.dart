// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';

void chkxer(
  final String SRNAMT,
  final int INFOT,
  final Nout NOUT,
  final Box<bool> LERR,
  final Box<bool> OK, [
  final TestDriver? test,
]) {
  final reason =
      'Illegal value of parameter number ${INFOT.i2} not detected by ${SRNAMT.trim().a6}';
  test?.expect(LERR.value, true, reason: reason);
  if (!LERR.value) {
    NOUT.println(' *** $reason ***');
    OK.value = false;
  }
  LERR.value = false;
}
