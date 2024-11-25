// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:test/test.dart';

import '../lin/alasum.dart';
import '../matgen/dlarnd.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dglmts.dart';
import 'dlatb9.dart';

Future<void> dckglm(
  final int NN,
  final Array<int> MVAL_,
  final Array<int> PVAL_,
  final Array<int> NVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AF_,
  final Array<double> B_,
  final Array<double> BF_,
  final Array<double> X_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
  final TestDriver test,
  final String group,
) async {
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final NVAL = NVAL_.having();
  final ISEED = ISEED_.having();
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final BF = BF_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTYPES = 8;
  final DOTYPE = Array<bool>(NTYPES);
  const PATH = 'GLM';

  // Initialize constants.

  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  final LDA = NMAX;
  final LDB = NMAX;
  final LWORK = NMAX * NMAX;

  // Check for valid input values.

  for (var IK = 1, FIRSTT = true; IK <= NN; IK++) {
    final M = MVAL[IK];
    final P = PVAL[IK];
    final N = NVAL[IK];
    if (M > N || N > M + P) {
      if (FIRSTT) {
        NOUT.println();
        FIRSTT = false;
      }
      NOUT.println(
          ' *** Invalid input  for GLM:  M = ${M.i6}, P = ${P.i6}, N = ${N.i6};\n'
          '     must satisfy M <= N <= M+P  (this set of values will be skipped)');
    }
  }

  test.group(group, () {
    var FIRSTT = true;

    // Do for each value of M in MVAL.
    for (final IK in 1.through(NN)) {
      final M = MVAL[IK];
      final P = PVAL[IK];
      final N = NVAL[IK];
      if (M > N || N > M + P) continue;

      for (final IMAT in 1.through(NTYPES)) {
        // Do the tests only if DOTYPE[ IMAT ] is true.
        final skip = !DOTYPE[IMAT];
        test('DCKGLM (M = $M, N = $N, P = $P TYPE = $IMAT)', () {
          final IINFO = Box(0);

          // Set up parameters with DLATB9 and generate test
          // matrices A and B with DLATMS.

          final (
            :TYPE,
            :KLA,
            :KUA,
            :KLB,
            :KUB,
            :ANORM,
            :BNORM,
            :MODEA,
            :MODEB,
            :CNDNMA,
            :CNDNMB,
            :DISTA,
            :DISTB
          ) = dlatb9(PATH, IMAT, M, P, N);

          dlatms(N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA,
              KUA, 'No packing', A.asMatrix(LDA), LDA, WORK, IINFO);
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9999(NOUT, IINFO.value);
            INFO.value = IINFO.value.abs();
            return;
          }

          dlatms(N, P, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB,
              KUB, 'No packing', B.asMatrix(LDB), LDB, WORK, IINFO);
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9999(NOUT, IINFO.value);
            INFO.value = IINFO.value.abs();
            return;
          }

          // Generate random left hand side vector of GLM

          for (var I = 1; I <= N; I++) {
            X[I] = dlarnd(2, ISEED);
          }

          final RESID = Box(0.0);
          dglmts(
              N,
              M,
              P,
              A.asMatrix(LDA),
              AF.asMatrix(LDA),
              LDA,
              B.asMatrix(LDB),
              BF.asMatrix(LDB),
              LDB,
              X,
              X(NMAX + 1),
              X(2 * NMAX + 1),
              X(3 * NMAX + 1),
              WORK,
              LWORK,
              RWORK,
              RESID);

          // Print information about the tests that did not
          // pass the threshold.

          final reason =
              ' N=${M.i4} M=${N.i4}, P=${P.i4}, type ${IMAT.i2}, test 1, ratio=${RESID.value.g13_6}';
          test.expect(RESID.value, lessThan(THRESH), reason: reason);
          if (RESID.value >= THRESH) {
            if (NFAIL == 0 && FIRSTT) {
              FIRSTT = false;
              alahdg(NOUT, PATH);
            }
            NOUT.println(reason);
            NFAIL++;
          }
          NRUN++;
        }, skip: skip);
      }
    }
  });

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGLM INFO = ${info.i5}');
}
