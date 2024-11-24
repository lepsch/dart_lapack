// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrlqt.dart';
import 'dlqt04.dart';

void dchklqt(
  final double THRESH,
  final bool TSTERR,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Nout NOUT,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  const NTESTS = 6;
  final RESULT = Array<double>(NTESTS);

  // Initialize constants

  final PATH = 'DTQ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);

  // Test the error exits
  test.group('error exits', () {
    if (TSTERR) derrlqt(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  // Do for each value of M in MVAL.
  for (final I in 1.through(NM)) {
    final M = MVAL[I];

    // Do for each value of N in NVAL.
    for (final J in 1.through(NN)) {
      final N = NVAL[J];

      // Do for each possible value of NB
      final MINMN = min(M, N);
      for (final K in 1.through(NNB)) {
        final NB = NBVAL[K];

        // Test DGELQT and DGEMLQT
        if ((NB <= MINMN) && (NB > 0)) {
          test('DGELQT and DGEMLQT (I=$I J=$J K=$K M=$M N=$N NB=$NB)', () {
            dlqt04(M, N, NB, RESULT);

            // Print information about the tests that did not
            // pass the threshold.

            for (var T = 1; T <= NTESTS; T++) {
              final reason =
                  ' M=${M.i5}, N=${N.i5}, NB=${NB.i4} test(${T.i2})=${RESULT[T].g12_5}';
              test.expect(RESULT[T], lessThan(THRESH), reason: reason);
              if (RESULT[T] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += NTESTS;
          });
        }
      }
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
