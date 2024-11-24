// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'zerrunhr_col.dart';
import 'zunhr_col01.dart';
import 'zunhr_col02.dart';

void zchkunhr_col(
  final double THRESH,
  final bool TSTERR,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Nout NOUT,
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

  final PATH = 'ZHH';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);

  // Test the error exits

  if (TSTERR) zerrunhr_col(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of M in MVAL.

  for (var I = 1; I <= NM; I++) {
    final M = MVAL[I];

    // Do for each value of N in NVAL.

    for (var J = 1; J <= NN; J++) {
      final N = NVAL[J];

      // Only for M >= N

      if (min(M, N) > 0 && M >= N) {
        // Do for each possible value of MB1

        for (var IMB1 = 1; IMB1 <= NNB; IMB1++) {
          final MB1 = NBVAL[IMB1];

          // Only for MB1 > N

          if (MB1 > N) {
            // Do for each possible value of NB1

            for (var INB1 = 1; INB1 <= NNB; INB1++) {
              final NB1 = NBVAL[INB1];

              // Do for each possible value of NB2

              for (var INB2 = 1; INB2 <= NNB; INB2++) {
                final NB2 = NBVAL[INB2];

                if (NB1 > 0 && NB2 > 0) {
                  // Test ZUNHR_COL

                  zunhr_col01(M, N, MB1, NB1, NB2, RESULT);

                  // Print information about the tests that did
                  // not pass the threshold.

                  for (var T = 1; T <= NTESTS; T++) {
                    if (RESULT[T] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                      NOUT.println(
                          'ZUNGTSQR and ZUNHR_COL: M=${M.i5}, N=${N.i5}, MB1=${MB1.i5}, NB1=${NB1.i5}, NB2=${NB2.i5} test(${T.i2})=${RESULT[T].g12_5}');
                      NFAIL++;
                    }
                  }
                  NRUN += NTESTS;
                }
              }
            }
          }
        }
      }
    }
  }

  // Do for each value of M in MVAL.

  for (var I = 1; I <= NM; I++) {
    final M = MVAL[I];

    // Do for each value of N in NVAL.

    for (var J = 1; J <= NN; J++) {
      final N = NVAL[J];

      // Only for M >= N

      if (min(M, N) > 0 && M >= N) {
        // Do for each possible value of MB1

        for (var IMB1 = 1; IMB1 <= NNB; IMB1++) {
          final MB1 = NBVAL[IMB1];

          // Only for MB1 > N

          if (MB1 > N) {
            // Do for each possible value of NB1

            for (var INB1 = 1; INB1 <= NNB; INB1++) {
              final NB1 = NBVAL[INB1];

              // Do for each possible value of NB2

              for (var INB2 = 1; INB2 <= NNB; INB2++) {
                final NB2 = NBVAL[INB2];

                if (NB1 > 0 && NB2 > 0) {
                  // Test ZUNHR_COL

                  zunhr_col02(M, N, MB1, NB1, NB2, RESULT);

                  // Print information about the tests that did
                  // not pass the threshold.

                  for (var T = 1; T <= NTESTS; T++) {
                    if (RESULT[T] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                      NOUT.println(
                          'ZUNGTSQR_ROW and ZUNHR_COL: M=${M.i5}, N=${N.i5}, MB1=${MB1.i5}, NB1=${NB1.i5}, NB2=${NB2.i5} test(${T.i2})=${RESULT[T].g12_5}');
                      NFAIL++;
                    }
                  }
                  NRUN += NTESTS;
                }
              }
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
