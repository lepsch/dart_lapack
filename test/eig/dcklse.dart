// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:test/test.dart';

import '../lin/alasum.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dlarhs.dart';
import 'dlatb9.dart';
import 'dlsets.dart';

Future<void> dcklse(
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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
  const NTESTS = 7;
  const NTYPES = 8;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  const PATH = 'LSE';

  // Initialize constants and the random number seed.

  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  final LDA = NMAX;
  final LDB = NMAX;
  final LWORK = NMAX * NMAX;

  test.group(group, () {
    var FIRSTT = true;

    // Check for valid input values.

    for (var IK = 1; IK <= NN; IK++) {
      final M = MVAL[IK];
      final P = PVAL[IK];
      final N = NVAL[IK];
      if (P > N || N > M + P) {
        if (FIRSTT) {
          NOUT.println();
          FIRSTT = false;
        }
        NOUT.println(
            ' *** Invalid input  for LSE:  M = ${M.i6}, P = ${P.i6}, N = ${N.i6};\n     must satisfy P <= N <= P+M  (this set of values will be skipped)');
      }
    }
    FIRSTT = true;

    // Do for each value of M in MVAL.

    for (final IK in 1.through(NN)) {
      final M = MVAL[IK];
      final P = PVAL[IK];
      final N = NVAL[IK];
      if (P > N || N > M + P) continue;

      for (final IMAT in 1.through(NTYPES)) {
        // Do the tests only if DOTYPE[ IMAT ] is true.
        final skip = !DOTYPE[IMAT];
        test('DCKLSE (M = $M, N = $N, P = $P TYPE = $IMAT)', () {
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

          dlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA,
              KUA, 'No packing', A.asMatrix(LDA), LDA, WORK, IINFO);
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9999(NOUT, IINFO.value);
            INFO.value = IINFO.value.abs();
            return;
          }

          dlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB,
              KUB, 'No packing', B.asMatrix(LDB), LDB, WORK, IINFO);
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9999(NOUT, IINFO.value);
            INFO.value = IINFO.value.abs();
            return;
          }

          // Generate the right-hand sides C and D for the LSE.
          var LD1 = max(N, 1), LD2 = max(M, 1);

          dlarhs(
              'DGE',
              'New solution',
              'Upper',
              'N',
              M,
              N,
              max(M - 1, 0),
              max(N - 1, 0),
              1,
              A.asMatrix(LDA),
              LDA,
              X(4 * NMAX + 1).asMatrix(LD1),
              LD1,
              X.asMatrix(LD2),
              LD2,
              ISEED,
              IINFO);

          LD1 = max(N, 1);
          LD2 = max(P, 1);
          dlarhs(
              'DGE',
              'Computed',
              'Upper',
              'N',
              P,
              N,
              max(P - 1, 0),
              max(N - 1, 0),
              1,
              B.asMatrix(LDB),
              LDB,
              X(4 * NMAX + 1).asMatrix(LD1),
              LD1,
              X(2 * NMAX + 1).asMatrix(LD2),
              LD2,
              ISEED,
              IINFO);

          const NT = 2;

          dlsets(
              M,
              P,
              N,
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
              X(4 * NMAX + 1),
              WORK,
              LWORK,
              RWORK,
              RESULT(1));

          // Print information about the tests that did not
          // pass the threshold.

          for (var I = 1; I <= NT; I++) {
            final reason =
                ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}';
            test.expect(RESULT[I], lessThan(THRESH), reason: reason);
            if (RESULT[I] >= THRESH) {
              if (NFAIL == 0 && FIRSTT) {
                FIRSTT = false;
                alahdg(NOUT, PATH);
              }
              NOUT.println(reason);
              NFAIL++;
            }
          }
          NRUN += NT;
        }, skip: skip);
      }
    }
  });

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout nout, final int info) {
  nout.println(' DLATMS in DCKLSE   INFO = ${info.i5}');
}
