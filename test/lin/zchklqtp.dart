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
import 'zerrlqtp.dart';
import 'zlqt05.dart';

void zchklqtp(
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
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  const NTESTS = 6;
  final RESULT = Array<double>(NTESTS);

  // Initialize constants

  const PATH = 'ZXQ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);

  // Test the error exits

  if (TSTERR) zerrlqtp(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of M

  for (var I = 1; I <= NM; I++) {
    final M = MVAL[I];

    // Do for each value of N

    for (var J = 1; J <= NN; J++) {
      final N = NVAL[J];

      // Do for each value of L

      final MINMN = min(M, N);
      for (var L = 0;
          max(MINMN, 1) < 0 ? L >= MINMN : L <= MINMN;
          L += max(MINMN, 1)) {
        // Do for each possible value of NB

        for (var K = 1; K <= NNB; K++) {
          final NB = NBVAL[K];

          // Test DTPLQT and DTPMLQT

          if ((NB <= M) && (NB > 0)) {
            zlqt05(M, N, L, NB, RESULT);

            // Print information about the tests that did not
            // pass the threshold.

            for (var T = 1; T <= NTESTS; T++) {
              if (RESULT[T] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' M=${M.i5}, N=${N.i5}, NB=${NB.i4} L=${L.i4} test(${T.i2})=${RESULT[T].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NTESTS;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
