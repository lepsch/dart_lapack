// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:test/test.dart';

import '../lin/alasum.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dgqrts.dart';
import 'dgrqts.dart';
import 'dlatb9.dart';

Future<void> dckgqr(
  final int NM,
  final Array<int> MVAL_,
  final int NP,
  final Array<int> PVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AF_,
  final Array<double> AQ_,
  final Array<double> AR_,
  final Array<double> TAUA_,
  final Array<double> B_,
  final Array<double> BF_,
  final Array<double> BZ_,
  final Array<double> BT_,
  final Array<double> BWK_,
  final Array<double> TAUB_,
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
  final AQ = AQ_.having();
  final AR = AR_.having();
  final TAUA = TAUA_.having();
  final B = B_.having();
  final BF = BF_.having();
  final BZ = BZ_.having();
  final BT = BT_.having();
  final BWK = BWK_.having();
  final TAUB = TAUB_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTESTS = 7;
  const NTYPES = 8;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  const PATH = 'GQR';

  // Initialize constants.

  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  final LDA = NMAX;
  final LDB = NMAX;
  final LWORK = NMAX * NMAX;

  test.group(group, () {
    var FIRSTT = true;

    // Do for each value of M in MVAL.
    for (final IM in 1.through(NM)) {
      final M = MVAL[IM];

      // Do for each value of P in PVAL.
      for (final IP in 1.through(NP)) {
        final P = PVAL[IP];

        // Do for each value of N in NVAL.
        for (final IN in 1.through(NN)) {
          final N = NVAL[IN];

          for (final IMAT in 1.through(NTYPES)) {
            // Do the tests only if DOTYPE[ IMAT ] is true.
            final skip = !DOTYPE[IMAT];
            test('DCKGQR (M = $M, N = $N, P = $P TYPE = $IMAT)', () {
              final IINFO = Box(0);

              // Test DGGRQF
              {
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
                ) = dlatb9('GRQ', IMAT, M, P, N);

                // Generate M by N matrix A

                dlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM,
                    KLA, KUA, 'No packing', A.asMatrix(LDA), LDA, WORK, IINFO);
                test.expect(IINFO.value, 0);
                if (IINFO.value != 0) {
                  print9999(NOUT, IINFO.value);
                  INFO.value = IINFO.value.abs();
                  return;
                }

                // Generate P by N matrix B

                dlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM,
                    KLB, KUB, 'No packing', B.asMatrix(LDB), LDB, WORK, IINFO);
                test.expect(IINFO.value, 0);
                if (IINFO.value != 0) {
                  print9999(NOUT, IINFO.value);
                  INFO.value = IINFO.value.abs();
                  return;
                }

                const NT = 4;

                dgrqts(
                    M,
                    P,
                    N,
                    A.asMatrix(LDA),
                    AF.asMatrix(LDA),
                    AQ.asMatrix(LDA),
                    AR.asMatrix(LDA),
                    LDA,
                    TAUA,
                    B.asMatrix(LDB),
                    BF.asMatrix(LDB),
                    BZ.asMatrix(LDB),
                    BT.asMatrix(LDB),
                    BWK.asMatrix(LDB),
                    LDB,
                    TAUB,
                    WORK,
                    LWORK,
                    RWORK,
                    RESULT);

                // Print information about the tests that did not
                // pass the threshold.

                for (var I = 1; I <= NT; I++) {
                  final reason =
                      ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}';
                  test.expect(RESULT[I], lessThan(THRESH), reason: reason);
                  if (RESULT[I] >= THRESH) {
                    if (NFAIL == 0 && FIRSTT) {
                      FIRSTT = false;
                      alahdg(NOUT, 'GRQ');
                    }
                    NOUT.println(reason);
                    NFAIL++;
                  }
                }
                NRUN += NT;
              }

              // Test DGGQRF
              {
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
                  MODEB: _,
                  :CNDNMA,
                  CNDNMB: _,
                  :DISTA,
                  :DISTB
                ) = dlatb9('GQR', IMAT, M, P, N);

                // Generate N-by-M matrix  A

                dlatms(N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM,
                    KLA, KUA, 'No packing', A.asMatrix(LDA), LDA, WORK, IINFO);
                test.expect(IINFO.value, 0);
                if (IINFO.value != 0) {
                  print9999(NOUT, IINFO.value);
                  INFO.value = IINFO.value.abs();
                  return;
                }

                // Generate N-by-P matrix  B

                dlatms(N, P, DISTB, ISEED, TYPE, RWORK, MODEA, CNDNMA, BNORM,
                    KLB, KUB, 'No packing', B.asMatrix(LDB), LDB, WORK, IINFO);
                test.expect(IINFO.value, 0);
                if (IINFO.value != 0) {
                  print9999(NOUT, IINFO.value);
                  INFO.value = IINFO.value.abs();
                  return;
                }

                const NT = 4;

                dgqrts(
                    N,
                    M,
                    P,
                    A.asMatrix(LDA),
                    AF.asMatrix(LDA),
                    AQ.asMatrix(LDA),
                    AR.asMatrix(LDA),
                    LDA,
                    TAUA,
                    B.asMatrix(LDB),
                    BF.asMatrix(LDB),
                    BZ.asMatrix(LDB),
                    BT.asMatrix(LDB),
                    BWK.asMatrix(LDB),
                    LDB,
                    TAUB,
                    WORK,
                    LWORK,
                    RWORK,
                    RESULT);

                // Print information about the tests that did not
                // pass the threshold.

                for (var I = 1; I <= NT; I++) {
                  final reason =
                      ' N=${N.i4} M=${M.i4}, P=${P.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}';
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
              }
            }, skip: skip);
          }
        }
      }
    }
  });

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGQR:    INFO = ${info.i5}');
}
