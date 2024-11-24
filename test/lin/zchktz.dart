// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgeqr2.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztzrzf.dart';

import '../matgen/zlatms.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dlaord.dart';
import 'zerrtz.dart';
import 'zqrt12.dart';
import 'zrzt01.dart';
import 'zrzt02.dart';

void zchktz(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final double THRESH,
  final bool TSTERR,
  final Array<Complex> A_,
  final Array<Complex> COPYA_,
  final Array<double> S_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
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
  final RWORK = RWORK_.having();

  const NTYPES = 3;
  const NTESTS = 3;
  const ONE = 1.0, ZERO = 0.0;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}TZ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  final EPS = dlamch('Epsilon');

  // Test the error exits

  if (TSTERR) zerrtz(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IM = 1; IM <= NM; IM++) {
    // Do for each value of M in MVAL.

    final M = MVAL[IM];
    final LDA = max(1, M);

    for (var IN = 1; IN <= NN; IN++) {
      // Do for each value of N in NVAL for which M <= N.

      final N = NVAL[IN];
      final MNMIN = min(M, N);
      final LWORK = max(1, N * N + 4 * M + N);

      if (M <= N) {
        for (var IMODE = 1; IMODE <= NTYPES; IMODE++) {
          if (!DOTYPE[IMODE]) continue;

          // Do for each type of singular value distribution.
          //    0:  zero matrix
          //    1:  one small singular value
          //    2:  exponential distribution

          final MODE = IMODE - 1;

          // Test ZTZRQF

          // Generate test matrix of size m by n using
          // singular value distribution indicated by `mode'.

          if (MODE == 0) {
            zlaset('Full', M, N, Complex.zero, Complex.zero, A.asMatrix(), LDA);
            for (var I = 1; I <= MNMIN; I++) {
              S[I] = ZERO;
            }
          } else {
            zlatms(M, N, 'Uniform', ISEED, 'Nonsymmetric', S, IMODE, ONE / EPS,
                ONE, M, N, 'No packing', A.asMatrix(), LDA, WORK, INFO);
            zgeqr2(M, N, A.asMatrix(), LDA, WORK, WORK(MNMIN + 1), INFO);
            zlaset('Lower', M - 1, N, Complex.zero, Complex.zero,
                A(2).asMatrix(), LDA);
            dlaord('Decreasing', MNMIN, S, 1);
          }

          // Save A and its singular values

          zlacpy('All', M, N, A.asMatrix(), LDA, COPYA.asMatrix(), LDA);

          // Call ZTZRZF to reduce the upper trapezoidal matrix to
          // upper triangular form.

          srnamc.SRNAMT = 'ZTZRZF';
          ztzrzf(M, N, A.asMatrix(), LDA, TAU, WORK, LWORK, INFO);

          // Compute norm(svd(a) - svd(r))

          RESULT[1] = zqrt12(M, M, A.asMatrix(), LDA, S, WORK, LWORK, RWORK);

          // Compute norm( A - R*Q )

          RESULT[2] = zrzt01(
              M, N, COPYA.asMatrix(), A.asMatrix(), LDA, TAU, WORK, LWORK);

          // Compute norm(Q'*Q - I).

          RESULT[3] = zrzt02(M, N, A.asMatrix(), LDA, TAU, WORK, LWORK);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NTESTS; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' M =${M.i5}, N =${N.i5}, type ${IMODE.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 3;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
