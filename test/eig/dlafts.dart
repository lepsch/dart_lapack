// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'dlahd2.dart';

void dlafts(
  final String TYPE,
  final int M,
  final int N,
  final int IMAT,
  final int NTESTS,
  final Array<double> RESULT_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout IOUNIT,
  final Box<int> IE, [
  final TestDriver? test,
]) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RESULT = RESULT_.having();
  final ISEED = ISEED_.having();

  if (M == N) {
    // Output for square matrices:

    for (var K = 1; K <= NTESTS; K++) {
      final reason = RESULT[K] < 10000.0
          ? ' Matrix order=${N.i5}, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')} result ${K.i3} is${RESULT[K].f8_2}'
          : ' Matrix order=${N.i5}, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')} result ${K.i3} is${(RESULT[K] * 10).d10_3}';

      test?.expect(RESULT[K], lessThan(THRESH), reason: reason);
      if (RESULT[K] >= THRESH) {
        // If this is the first test to fail, call DLAHD2
        // to print a header to the data file.
        if (IE.value == 0) dlahd2(IOUNIT, TYPE);
        IE.value++;
        IOUNIT.println(reason);
      }
    }
  } else {
    // Output for rectangular matrices

    for (var K = 1; K <= NTESTS; K++) {
      final reason = RESULT[K] < 10000.0
          ? ' ${M.i5} x${N.i5} matrix, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')}: result ${K.i3} is${RESULT[K].f8_2}'
          : ' ${M.i5} x${N.i5} matrix, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')}: result ${K.i3} is${(RESULT[K] * 10).d10_3}';

      test?.expect(RESULT[K], lessThan(THRESH), reason: reason);
      if (RESULT[K] >= THRESH) {
        // If this is the first test to fail, call DLAHD2
        // to print a header to the data file.
        if (IE.value == 0) dlahd2(IOUNIT, TYPE);
        IE.value++;
        IOUNIT.println(reason);
      }
    }
  }
}
