// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrtz.dart';
import 'dlaord.dart';
import 'dqrt12.dart';
import 'drzt01.dart';
import 'drzt02.dart';

void dchktz(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final double THRESH,
  final bool TSTERR,
  final Array<double> A_,
  final Array<double> COPYA_,
  final Array<double> S_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Nout NOUT,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final COPYA = COPYA_.having();
  final S = S_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();

  const NTYPES = 3, NTESTS = 3;
  const ONE = 1.0, ZERO = 0.0;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}TZ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);
  final EPS = dlamch('Epsilon');

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrtz(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  for (final IM in 1.through(NM)) {
    // Do for each value of M in MVAL.
    final M = MVAL[IM];
    final LDA = max(1, M).toInt();

    for (final IN in 1.through(NN)) {
      // Do for each value of N in NVAL for which M <= N.
      final N = NVAL[IN];
      final MNMIN = min(M, N);
      final LWORK = max(1, max(N * N + 4 * M + N, M * N + 2 * MNMIN + 4 * N));
      if (M > N) continue;

      for (final IMODE in 1.through(NTYPES)) {
        final skip = !DOTYPE[IMODE];
        test('DCHKTZ (IM=$IM IN=$IN IMODE=$IMODE)', () {
          final INFO = Box(0);

          // Do for each type of singular value distribution.
          //    0:  zero matrix
          //    1:  one small singular value
          //    2:  exponential distribution
          final MODE = IMODE - 1;

          // Test DTZRQF

          // Generate test matrix of size m by n using
          // singular value distribution indicated by `mode'.
          if (MODE == 0) {
            dlaset('Full', M, N, ZERO, ZERO, A.asMatrix(), LDA);
            for (var I = 1; I <= MNMIN; I++) {
              S[I] = ZERO;
            }
          } else {
            dlatms(M, N, 'Uniform', ISEED, 'Nonsymmetric', S, IMODE, ONE / EPS,
                ONE, M, N, 'No packing', A.asMatrix(), LDA, WORK, INFO);
            dgeqr2(M, N, A.asMatrix(), LDA, WORK, WORK(MNMIN + 1), INFO);
            dlaset('Lower', M - 1, N, ZERO, ZERO, A(2).asMatrix(), LDA);
            dlaord('Decreasing', MNMIN, S, 1);
          }

          // Save A and its singular values
          dlacpy('All', M, N, A.asMatrix(), LDA, COPYA.asMatrix(), LDA);

          // Call DTZRZF to reduce the upper trapezoidal matrix to
          // upper triangular form.
          srnamc.SRNAMT = 'DTZRZF';
          dtzrzf(M, N, A.asMatrix(), LDA, TAU, WORK, LWORK, INFO);

          // Compute norm(svd(a) - svd(r))
          RESULT[1] = dqrt12(M, M, A.asMatrix(), LDA, S, WORK, LWORK);

          // Compute norm( A - R*Q )
          RESULT[2] = drzt01(
              M, N, COPYA.asMatrix(), A.asMatrix(), LDA, TAU, WORK, LWORK);

          // Compute norm(Q'*Q - I).
          RESULT[3] = drzt02(M, N, A.asMatrix(), LDA, TAU, WORK, LWORK);

          // Print information about the tests that did not pass the threshold.
          for (var K = 1; K <= NTESTS; K++) {
            final reason =
                ' M =${M.i5}, N =${N.i5}, type ${IMODE.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}';
            test.expect(RESULT[K], lessThan(THRESH), reason: reason);
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(reason);
              NFAIL++;
            }
          }
          NRUN += 3;
        }, skip: skip);
      }
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
