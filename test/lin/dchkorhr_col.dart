// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrorhr_col.dart';
import 'dorhr_col01.dart';
import 'dorhr_col02.dart';

void dchkorhr_col(
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
  final PATH = 'DHH';
  var NRUN = 0;
  var NFAIL = 0;
  const NERRS = 0;

  // Test the error exits
  test.group('error exits', () {
    if (TSTERR) derrorhr_col(PATH, NOUT, test);
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

      // Only for M >= N
      if (min(M, N) <= 0 || M < N) continue;

      // Do for each possible value of MB1
      for (final IMB1 in 1.through(NNB)) {
        final MB1 = NBVAL[IMB1];

        // Only for MB1 > N
        if (MB1 <= N) continue;

        // Do for each possible value of NB1
        for (final INB1 in 1.through(NNB)) {
          final NB1 = NBVAL[INB1];

          // Do for each possible value of NB2
          for (final INB2 in 1.through(NNB)) {
            final NB2 = NBVAL[INB2];

            if (NB1 <= 0 || NB2 <= 0) continue;

            // Test DORHR_COL
            test(
                'DORHR_COL (I=$I J=$J IMB1=$IMB1 INB1=$INB1 INB2=$INB2 M=$M N=$N MB1=$MB1 NB1=$NB1 NB2=$NB2)',
                () {
              dorhr_col01(M, N, MB1, NB1, NB2, RESULT);

              // Print information about the tests that did
              // not pass the threshold.
              for (var T = 1; T <= NTESTS; T++) {
                final reason =
                    'DORGTSQR and DORHR_COL: M=${M.i5}, N=${N.i5}, MB1=${MB1.i5}, NB1=${NB1.i5}, NB2=${NB2.i5} test(${T.i2})=${RESULT[T].g12_5}';
                test.expect(RESULT[T], lessThan(THRESH), reason: reason);
                if (RESULT[T] >= THRESH) {
                  if (NFAIL == 0 && NERRS == 0) alahd(NOUT, PATH);
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
  }

  // Do for each value of M in MVAL.
  for (final I in 1.through(NM)) {
    final M = MVAL[I];

    // Do for each value of N in NVAL.
    for (final J in 1.through(NN)) {
      final N = NVAL[J];

      // Only for M >= N
      if (min(M, N) <= 0 || M < N) continue;

      // Do for each possible value of MB1
      for (final IMB1 in 1.through(NNB)) {
        final MB1 = NBVAL[IMB1];

        // Only for MB1 > N
        if (MB1 <= N) continue;

        // Do for each possible value of NB1
        for (final INB1 in 1.through(NNB)) {
          final NB1 = NBVAL[INB1];

          // Do for each possible value of NB2
          for (final INB2 in 1.through(NNB)) {
            final NB2 = NBVAL[INB2];

            if (NB1 <= 0 || NB2 <= 0) continue;

            // Test DORHR_COL
            test(
                'DORHR_COL02 (I=$I J=$J IMB1=$IMB1 INB1=$INB1 INB2=$INB2 M=$M N=$N MB1=$MB1 NB1=$NB1 NB2=$NB2)',
                () {
              dorhr_col02(M, N, MB1, NB1, NB2, RESULT);

              // Print information about the tests that did
              // not pass the threshold.
              for (var T = 1; T <= NTESTS; T++) {
                final reason =
                    'DORGTSQR_ROW and DORHR_COL: M=${M.i5}, N=${N.i5}, MB1=${MB1.i5}, NB1=${NB1.i5}, NB2=${NB2.i5} test(${T.i2})=${RESULT[T].g12_5}';
                test.expect(RESULT[T], lessThan(THRESH), reason: reason);
                if (RESULT[T] >= THRESH) {
                  if (NFAIL == 0 && NERRS == 0) alahd(NOUT, PATH);
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
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS);
}
