import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zerrtsqr.dart';
import 'ztsqr01.dart';

void zchktsqr(
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

  final PATH = 'ZTS';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);

  // Test the error exits

  xlaenv(1, 0);
  xlaenv(2, 0);
  if (TSTERR) zerrtsqr(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of M in MVAL.

  for (var I = 1; I <= NM; I++) {
    final M = MVAL[I];

    // Do for each value of N in NVAL.

    for (var J = 1; J <= NN; J++) {
      final N = NVAL[J];
      if (min(M, N) != 0) {
        for (var INB = 1; INB <= NNB; INB++) {
          final MB = NBVAL[INB];
          xlaenv(1, MB);
          for (var IMB = 1; IMB <= NNB; IMB++) {
            final NB = NBVAL[IMB];
            xlaenv(2, NB);

            // Test ZGEQR and ZGEMQR

            ztsqr01('TS', M, N, MB, NB, RESULT);

            // Print information about the tests that did not
            // pass the threshold.

            for (var T = 1; T <= NTESTS; T++) {
              if (RESULT[T] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    'TS: M=${M.i5}, N=${N.i5}, MB=${MB.i5}, NB=${NB.i5} test(${T.i2})=${RESULT[T].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NTESTS;
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
      if (min(M, N) != 0) {
        for (var INB = 1; INB <= NNB; INB++) {
          final MB = NBVAL[INB];
          xlaenv(1, MB);
          for (var IMB = 1; IMB <= NNB; IMB++) {
            final NB = NBVAL[IMB];
            xlaenv(2, NB);

            // Test ZGELQ and ZGEMLQ

            ztsqr01('SW', M, N, MB, NB, RESULT);

            // Print information about the tests that did not
            // pass the threshold.

            for (var T = 1; T <= NTESTS; T++) {
              if (RESULT[T] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    'SW: M=${M.i5}, N=${N.i5}, MB=${MB.i5}, NB=${NB.i5} test(${T.i2})=${RESULT[T].g12_5}');
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
