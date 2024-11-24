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
import 'package:lapack/src/zgeqp3.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatms.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dlaord.dart';
import 'icopy.dart';
import 'xlaenv.dart';
import 'zqpt01.dart';
import 'zqrt11.dart';
import 'zqrt12.dart';

void zchkq3(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Array<int> NXVAL_,
  final double THRESH,
  final Array<Complex> A_,
  final Array<Complex> COPYA_,
  final Array<double> S_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final NXVAL = NXVAL_.having();
  final A = A_.having();
  final COPYA = COPYA_.having();
  final S = S_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const NTYPES = 6, NTESTS = 3;
  const ONE = 1.0, ZERO = 0.0;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}Q3';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  final EPS = dlamch('Epsilon');
  infoc.INFOT = 0;

  for (var IM = 1; IM <= NM; IM++) {
    // Do for each value of M in MVAL.

    final M = MVAL[IM];
    final LDA = max(1, M);

    for (var IN = 1; IN <= NN; IN++) {
      // Do for each value of N in NVAL.

      final N = NVAL[IN];
      final MNMIN = min(M, N);
      final LWORK = max(1, (M * max(M, N) + 4 * MNMIN + max(M, N)).toInt());

      for (var IMODE = 1; IMODE <= NTYPES; IMODE++) {
        if (!DOTYPE[IMODE]) continue;

        // Do for each type of matrix
        //    1:  zero matrix
        //    2:  one small singular value
        //    3:  geometric distribution of singular values
        //    4:  first n/2 columns fixed
        //    5:  last n/2 columns fixed
        //    6:  every second column fixed

        final MODE = IMODE > 3 ? 1 : IMODE;

        // Generate test matrix of size m by n using
        // singular value distribution indicated by `mode'.

        for (var I = 1; I <= N; I++) {
          IWORK[I] = 0;
        }
        if (IMODE == 1) {
          zlaset(
              'Full', M, N, Complex.zero, Complex.zero, COPYA.asMatrix(), LDA);
          for (var I = 1; I <= MNMIN; I++) {
            S[I] = ZERO;
          }
        } else {
          zlatms(M, N, 'Uniform', ISEED, 'Nonsymm', S, MODE, ONE / EPS, ONE, M,
              N, 'No packing', COPYA.asMatrix(), LDA, WORK, INFO);
          if (IMODE >= 4) {
            final (ILOW, ISTEP, IHIGH) = switch (IMODE) {
              4 => (1, 1, max(1, N ~/ 2)),
              5 => (max(1, N ~/ 2), 1, N),
              6 => (1, 2, N),
              _ => throw UnimplementedError(),
            };
            for (var I = ILOW; I <= IHIGH; I += ISTEP) {
              IWORK[I] = 1;
            }
          }
          dlaord('Decreasing', MNMIN, S, 1);
        }

        for (var INB = 1; INB <= NNB; INB++) {
          // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

          final NB = NBVAL[INB];
          xlaenv(1, NB);
          final NX = NXVAL[INB];
          xlaenv(3, NX);

          // Save A and its singular values and a copy of
          // vector IWORK.

          zlacpy('All', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
          icopy(N, IWORK(1), 1, IWORK(N + 1), 1);

          // Workspace needed.

          final LW = NB * (N + 1);

          srnamc.SRNAMT = 'ZGEQP3';
          zgeqp3(M, N, A.asMatrix(), LDA, IWORK(N + 1), TAU, WORK, LW, RWORK,
              INFO);

          // Compute norm(svd(a) - svd(r))

          RESULT[1] = zqrt12(M, N, A.asMatrix(), LDA, S, WORK, LWORK, RWORK);

          // Compute norm( A*P - Q*R )

          RESULT[2] = zqpt01(M, N, MNMIN, COPYA.asMatrix(), A.asMatrix(), LDA,
              TAU, IWORK(N + 1), WORK, LWORK);

          // Compute Q'*Q

          RESULT[3] = zqrt11(M, MNMIN, A.asMatrix(), LDA, TAU, WORK, LWORK);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NTESTS; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' ZGEQP3 M =${M.i5}, N =${N.i5}, NB =${NB.i4}, type ${IMODE.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += NTESTS;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
