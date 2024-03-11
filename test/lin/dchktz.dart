import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqr2.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtzrzf.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
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
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}TZ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  final EPS = dlamch('Epsilon');

  // Test the error exits

  if (TSTERR) derrtz(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IM = 1; IM <= NM; IM++) {
    // Do for each value of M in MVAL.

    final M = MVAL[IM];
    final LDA = max(1, M).toInt();

    for (var IN = 1; IN <= NN; IN++) {
      // Do for each value of N in NVAL for which M <= N.

      final N = NVAL[IN];
      final MNMIN = min(M, N);
      final LWORK = max(1, max(N * N + 4 * M + N, M * N + 2 * MNMIN + 4 * N));

      if (M <= N) {
        for (var IMODE = 1; IMODE <= NTYPES; IMODE++) {
          if (!DOTYPE[IMODE]) continue;

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

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NTESTS; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' M =${M.i5}, N =${N.i5}, type ${IMODE.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL = NFAIL + 1;
            }
          }
          NRUN = NRUN + 3;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
