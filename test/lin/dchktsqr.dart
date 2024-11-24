// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrtsqr.dart';
import 'dtsqr01.dart';
import 'xlaenv.dart';

void dchktsqr(
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
  const PATH = 'DTS';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = 0;

  // Test the error exits
  test.group('error exits', () {
    test.setUp(() {
      xlaenv(1, 0);
      xlaenv(2, 0);
    });
    if (TSTERR) derrtsqr(PATH, NOUT, test);
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
      if (min(M, N) == 0) continue;

      for (final INB in 1.through(NNB)) {
        final MB = NBVAL[INB];

        for (final IMB in 1.through(NNB)) {
          final NB = NBVAL[IMB];

          // Test DGEQR and DGEMQR
          test(
              'DGEQR and DGEMQR (I=$I J=$J INB=$INB IMB=$IMB M=$M N=$N MB=$MB NB=$NB)',
              () {
            xlaenv(1, MB);
            xlaenv(2, NB);

            dtsqr01('TS', M, N, MB, NB, RESULT);

            // Print information about the tests that did not
            // pass the threshold.
            for (var T = 1; T <= NTESTS; T++) {
              final reason =
                  'TS: M=${M.i5}, N=${N.i5}, MB=${MB.i5}, NB=${NB.i5} test(${T.i2})=${RESULT[T].g12_5}';
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

  // Do for each value of M in MVAL.
  for (final I in 1.through(NM)) {
    final M = MVAL[I];

    // Do for each value of N in NVAL.
    for (final J in 1.through(NN)) {
      final N = NVAL[J];
      if (min(M, N) == 0) continue;

      for (final INB in 1.through(NNB)) {
        final MB = NBVAL[INB];

        for (final IMB in 1.through(NNB)) {
          final NB = NBVAL[IMB];

          // Test DGEQR and DGEMQR
          test(
              'DGEQR and DGEMQR (I=$I J=$J INB=$INB IMB=$IMB M=$M N=$N MB=$MB NB=$NB)',
              () {
            xlaenv(1, MB);
            xlaenv(2, NB);

            dtsqr01('SW', M, N, MB, NB, RESULT);

            // Print information about the tests that did not
            // pass the threshold.
            for (var T = 1; T <= NTESTS; T++) {
              final reason =
                  'SW: M=${M.i5}, N=${N.i5}, MB=${MB.i5}, NB=${NB.i5} test(${T.i2})=${RESULT[T].g12_5}';
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

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS);
}
