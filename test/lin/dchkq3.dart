import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqp3.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dlaord.dart';
import 'dqpt01.dart';
import 'dqrt11.dart';
import 'dqrt12.dart';
import 'icopy.dart';
import 'xlaenv.dart';

void dchkq3(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Array<int> NXVAL_,
  final double THRESH,
  final Array<double> A_,
  final Array<double> COPYA_,
  final Array<double> S_,
  final Array<double> TAU_,
  final Array<double> WORK_,
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
  final IWORK = IWORK_.having();
  const NTYPES = 6;
  const NTESTS = 3;
  const ONE = 1.0, ZERO = 0.0;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}Q3';
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
      final LWORK = max(
              1,
              max(M * max(M, N) + 4 * MNMIN + max(M, N),
                  M * N + 2 * MNMIN + 4 * N))
          .toInt();

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
          dlaset('Full', M, N, ZERO, ZERO, COPYA.asMatrix(), LDA);
          for (var I = 1; I <= MNMIN; I++) {
            S[I] = ZERO;
          }
        } else {
          dlatms(M, N, 'Uniform', ISEED, 'Nonsymm', S, MODE, ONE / EPS, ONE, M,
              N, 'No packing', COPYA.asMatrix(), LDA, WORK, INFO);
          if (IMODE >= 4) {
            var ILOW = 0, ISTEP = 0, IHIGH = 0;
            if (IMODE == 4) {
              ILOW = 1;
              ISTEP = 1;
              IHIGH = max(1, N ~/ 2);
            } else if (IMODE == 5) {
              ILOW = max(1, N ~/ 2);
              ISTEP = 1;
              IHIGH = N;
            } else if (IMODE == 6) {
              ILOW = 1;
              ISTEP = 2;
              IHIGH = N;
            }
            for (var I = ILOW;
                ISTEP < 0 ? I >= IHIGH : I <= IHIGH;
                I += ISTEP) {
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

          // Get a working copy of COPYA into A and a copy of
          // vector IWORK.

          dlacpy('All', M, N, COPYA.asMatrix(), LDA, A.asMatrix(), LDA);
          icopy(N, IWORK(1), 1, IWORK(N + 1), 1);

          // Compute the QR factorization with pivoting of A

          final LW = max(1, 2 * N + NB * (N + 1));

          // Compute the QP3 factorization of A

          srnamc.SRNAMT = 'DGEQP3';
          dgeqp3(M, N, A.asMatrix(), LDA, IWORK(N + 1), TAU, WORK, LW, INFO);

          // Compute norm(svd(a) - svd(r))

          RESULT[1] = dqrt12(M, N, A.asMatrix(), LDA, S, WORK, LWORK);

          // Compute norm( A*P - Q*R )

          RESULT[2] = dqpt01(M, N, MNMIN, COPYA.asMatrix(), A.asMatrix(), LDA,
              TAU, IWORK(N + 1), WORK, LWORK);

          // Compute Q'*Q

          RESULT[3] = dqrt11(M, MNMIN, A.asMatrix(), LDA, TAU, WORK, LWORK);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NTESTS; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' DGEQP3 M =${M.i5}, N =${N.i5}, NB =${NB.i4}, type ${IMODE.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
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
