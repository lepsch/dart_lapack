// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dsgesv.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'common.dart';
import 'dget08.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';

void ddrvab(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Array<double> SWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final B = B_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final SWORK = SWORK_.having();
  final IWORK = IWORK_.having();

  const ZERO = 0.0;
  const NTYPES = 11;
  const NTESTS = 1;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);

  const ISEEDY = [2006, 2007, 2008, 2009];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}GE';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  infoc.INFOT = 0;

  // Do for each value of M in MVAL

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);

    final N = M;
    final NIMAT = M <= 0 || N <= 0 ? 1 : NTYPES;

    imatLoop:
    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 5, 6, or 7 if the matrix size is too small.

      final ZEROT = IMAT >= 5 && IMAT <= 7;
      if (ZEROT && N < IMAT - 4) continue;

      // Set up parameters with DLATB4 and generate a test matrix
      // with DLATMS.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
          dlatb4(PATH, IMAT, M, N);

      srnamc.SRNAMT = 'DLATMS';
      dlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
          'No packing', A.asMatrix(), LDA, WORK, INFO);

      // Check error code from DLATMS.

      if (INFO.value != 0) {
        alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
        continue;
      }

      // For types 5-7, zero one or more columns of the matrix to
      // test that INFO is returned correctly.

      final int IZERO;
      if (ZEROT) {
        if (IMAT == 5) {
          IZERO = 1;
        } else if (IMAT == 6) {
          IZERO = min(M, N);
        } else {
          IZERO = min(M, N) ~/ 2 + 1;
        }
        var IOFF = (IZERO - 1) * LDA;
        if (IMAT < 7) {
          for (var I = 1; I <= M; I++) {
            A[IOFF + I] = ZERO;
          }
        } else {
          dlaset('Full', M, N - IZERO + 1, ZERO, ZERO, A(IOFF + 1).asMatrix(),
              LDA);
        }
      } else {
        IZERO = 0;
      }

      for (var IRHS = 1; IRHS <= NNS; IRHS++) {
        final NRHS = NSVAL[IRHS];
        final XTYPE = 'N';
        final TRANS = 'N';

        srnamc.SRNAMT = 'DLARHS';
        dlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A.asMatrix(), LDA,
            X.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);

        srnamc.SRNAMT = 'DSGESV';

        dlacpy('Full', M, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

        final ITER = Box(0);
        dsgesv(N, NRHS, A.asMatrix(), LDA, IWORK, B.asMatrix(), LDA,
            X.asMatrix(), LDA, WORK.asMatrix(), SWORK, ITER, INFO);

        if (ITER.value < 0) {
          dlacpy('Full', M, N, AFAC.asMatrix(), LDA, A.asMatrix(), LDA);
        }

        // Check error code from DSGESV. This should be the same as
        // the one of DGETRF.

        if (INFO.value != IZERO) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NERRS.value++;

          if (INFO.value != IZERO && IZERO != 0) {
            NOUT.println(
                ' *** DSGESV returned with INFO =${INFO.value.i5} instead of ${IZERO.i5}\n ==> M =${M.i5}, type ${IMAT.i2}');
          } else {
            NOUT.println(
                ' *** Error code from DSGESV=${INFO.value.i5} for M=${M.i5}, type ${IMAT.i2}');
          }
        }

        // Skip the remaining test if the matrix is singular.

        if (INFO.value != 0) continue imatLoop;

        // Check the quality of the solution

        dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

        dget08(TRANS, N, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
            WORK.asMatrix(), LDA, RWORK, RESULT(1));

        // Check if the test passes the testing.
        // Print information about the tests that did not
        // pass the testing.

        // If iterative refinement has been used and claimed to
        // be successful (ITER>0), we want
        //   NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1

        // If double precision has been used (ITER<0), we want
        //   NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
        // (Cf. the linear solver testing routines)

        if ((THRESH <= 0.0e+00) ||
            ((ITER.value >= 0) && (N > 0) && (RESULT[1] >= sqrt(N))) ||
            ((ITER.value < 0) && (RESULT[1] >= THRESH))) {
          if (NFAIL == 0 && NERRS.value == 0) {
            NOUT.println('\n DGE:  General dense matrices');
            NOUT.println(' Matrix types:');
            NOUT.println(
                '    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero');
            NOUT.println(' Test ratios:');
            NOUT.println(
                '   ${1.i2}: norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF\n    or norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF');
            NOUT.println(' Messages:');
          }

          NOUT.println(
              ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${1.i2}) =${RESULT[1].g12_5}');
          NFAIL++;
        }
        NRUN++;
      }
    }
  }

  // Print a summary of the results.

  if (NFAIL > 0) {
    NOUT.println(
        ' DSGESV: ${NFAIL.i6} out of ${NRUN.i6} tests failed to pass the threshold');
  } else {
    NOUT.println(
        '\n All tests for DSGESV routines passed the threshold ( ${NRUN.i6} tests run)');
  }
  if (NERRS.value > 0) {
    NOUT.println('${' ' * 6}${NERRS.value.i6} error messages recorded');
  }
}
