// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrtr.dart';
import 'dget04.dart';
import 'dlarhs.dart';
import 'dlattb.dart';
import 'dtbt02.dart';
import 'dtbt03.dart';
import 'dtbt05.dart';
import 'dtbt06.dart';

void dchktb(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> AB_,
  final Array<double> AINV_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
  final TestDriver test,
) {
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NSVAL = NSVAL_.having();
  final AB = AB_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const NTYPE1 = 9, NTYPES = 17;
  const NTESTS = 8;
  const NTRAN = 3;
  const ONE = 1.0, ZERO = 0.0;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], TRANSS = ['N', 'T', 'C'];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}TB';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrtr(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  // Do for each value of N in NVAL
  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPE1;
    final NIMAT2 = N <= 0 ? NTYPE1 + 1 : NTYPES;

    final NK = min(N + 1, 4);
    for (final IK in 1.through(NK)) {
      // Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
      // it easier to skip redundant values for small values of N.

      final KD = switch (IK) {
        2 => max(N, 0),
        3 => (3 * N - 1) ~/ 4,
        4 => (N + 1) ~/ 4,
        1 || _ => 0,
      };
      final LDAB = KD + 1;

      for (final IMAT in 1.through(NIMAT)) {
        // Do the tests only if DOTYPE( IMAT ) is true.
        final skip = !DOTYPE[IMAT];

        test('DCHKTB - 1 (IN=$IN NK=$NK IMAGE=$IMAT)', () {
          final INFO = Box(0);

          for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
            // Do first for UPLO = 'U', then for UPLO = 'L'
            final UPLO = UPLOS[IUPLO - 1];

            // Call DLATTB to generate a triangular test matrix.
            final DIAG = Box('');
            srnamc.SRNAMT = 'DLATTB';
            dlattb(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, KD,
                AB.asMatrix(), LDAB, X, WORK, INFO);

            // Set IDIAG = 1 for non-unit matrices, 2 for unit.
            final IDIAG = lsame(DIAG.value, 'N') ? 1 : 2;

            // Form the inverse of A so we can get a good estimate
            // of RCONDC = 1/(norm(A) * norm(inv(A))).
            dlaset('Full', N, N, ZERO, ONE, AINV.asMatrix(), LDA);
            if (lsame(UPLO, 'U')) {
              for (var J = 1; J <= N; J++) {
                dtbsv(UPLO, 'No transpose', DIAG.value, J, KD, AB.asMatrix(),
                    LDAB, AINV((J - 1) * LDA + 1), 1);
              }
            } else {
              for (var J = 1; J <= N; J++) {
                dtbsv(
                    UPLO,
                    'No transpose',
                    DIAG.value,
                    N - J + 1,
                    KD,
                    AB((J - 1) * LDAB + 1).asMatrix(),
                    LDAB,
                    AINV((J - 1) * LDA + J),
                    1);
              }
            }

            // Compute the 1-norm condition number of A.
            var ANORM = dlantb(
                '1', UPLO, DIAG.value, N, KD, AB.asMatrix(), LDAB, RWORK);
            var AINVNM = dlantr(
                '1', UPLO, DIAG.value, N, N, AINV.asMatrix(), LDA, RWORK);
            final double RCONDO;
            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDO = ONE;
            } else {
              RCONDO = (ONE / ANORM) / AINVNM;
            }

            // Compute the infinity-norm condition number of A.
            ANORM = dlantb(
                'I', UPLO, DIAG.value, N, KD, AB.asMatrix(), LDAB, RWORK);
            AINVNM = dlantr(
                'I', UPLO, DIAG.value, N, N, AINV.asMatrix(), LDA, RWORK);
            final double RCONDI;
            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDI = ONE;
            } else {
              RCONDI = (ONE / ANORM) / AINVNM;
            }

            for (var IRHS = 1; IRHS <= NNS; IRHS++) {
              final NRHS = NSVAL[IRHS];
              var XTYPE = 'N';

              for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
                // Do for op(A) = A, A**T, or A**H.

                final TRANS = TRANSS[ITRAN - 1];
                final (_, RCONDC) = ITRAN == 1 ? ('O', RCONDO) : ('I', RCONDI);

                // +    TEST 1
                // Solve and compute residual for op(A)*x = b.

                srnamc.SRNAMT = 'DLARHS';
                dlarhs(
                    PATH,
                    XTYPE,
                    UPLO,
                    TRANS,
                    N,
                    N,
                    KD,
                    IDIAG,
                    NRHS,
                    AB.asMatrix(),
                    LDAB,
                    XACT.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDA,
                    ISEED,
                    INFO);
                XTYPE = 'C';
                dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                srnamc.SRNAMT = 'DTBTRS';
                dtbtrs(UPLO, TRANS, DIAG.value, N, KD, NRHS, AB.asMatrix(),
                    LDAB, X.asMatrix(), LDA, INFO);

                // Check error code from DTBTRS.
                test.expect(INFO.value, 0);
                if (INFO.value != 0) {
                  alaerh(
                      PATH,
                      'DTBTRS',
                      INFO.value,
                      0,
                      UPLO + TRANS + DIAG.value,
                      N,
                      N,
                      KD,
                      KD,
                      NRHS,
                      IMAT,
                      NFAIL,
                      NERRS,
                      NOUT);
                }

                dtbt02(
                    UPLO,
                    TRANS,
                    DIAG.value,
                    N,
                    KD,
                    NRHS,
                    AB.asMatrix(),
                    LDAB,
                    X.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDA,
                    WORK,
                    RESULT(1));

                // +    TEST 2
                // Check solution from generated exact solution.

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(2));

                // +    TESTS 3, 4, and 5
                // Use iterative refinement to improve the solution
                // and compute error bounds.

                srnamc.SRNAMT = 'DTBRFS';
                dtbrfs(
                    UPLO,
                    TRANS,
                    DIAG.value,
                    N,
                    KD,
                    NRHS,
                    AB.asMatrix(),
                    LDAB,
                    B.asMatrix(),
                    LDA,
                    X.asMatrix(),
                    LDA,
                    RWORK,
                    RWORK(NRHS + 1),
                    WORK,
                    IWORK,
                    INFO);

                // Check error code from DTBRFS.
                test.expect(INFO.value, 0);
                if (INFO.value != 0) {
                  alaerh(
                      PATH,
                      'DTBRFS',
                      INFO.value,
                      0,
                      UPLO + TRANS + DIAG.value,
                      N,
                      N,
                      KD,
                      KD,
                      NRHS,
                      IMAT,
                      NFAIL,
                      NERRS,
                      NOUT);
                }

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                dtbt05(
                    UPLO,
                    TRANS,
                    DIAG.value,
                    N,
                    KD,
                    NRHS,
                    AB.asMatrix(),
                    LDAB,
                    B.asMatrix(),
                    LDA,
                    X.asMatrix(),
                    LDA,
                    XACT.asMatrix(),
                    LDA,
                    RWORK,
                    RWORK(NRHS + 1),
                    RESULT(4));

                // Print information about the tests that did not
                // pass the threshold.

                for (var K = 1; K <= 5; K++) {
                  final reason =
                      ' UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.value.a1}\', N=${N.i5}, KD=${KD.i5}, NRHS=${NRHS.i5}, type ${IMAT.i2}, test(${K.i2})=${RESULT[K].g12_5}';
                  test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.println(reason);
                    NFAIL++;
                  }
                }
                NRUN += 5;
              }
            }

            // +    TEST 6
            // Get an estimate of RCOND = 1/CNDNUM.

            for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
              final (NORM, RCONDC) = ITRAN == 1 ? ('O', RCONDO) : ('I', RCONDI);
              final RCOND = Box(ZERO);
              srnamc.SRNAMT = 'DTBCON';
              dtbcon(NORM, UPLO, DIAG.value, N, KD, AB.asMatrix(), LDAB, RCOND,
                  WORK, IWORK, INFO);

              // Check error code from DTBCON.
              test.expect(INFO.value, 0);
              if (INFO.value != 0) {
                alaerh(PATH, 'DTBCON', INFO.value, 0, NORM + UPLO + DIAG.value,
                    N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT);
              }

              dtbt06(RCOND.value, RCONDC, UPLO, DIAG.value, N, KD,
                  AB.asMatrix(), LDAB, RWORK, RESULT(6));

              // Print information about the tests that did not pass
              // the threshold.
              final reason =
                  ' DTBCON( \'${NORM.a1}\'${UPLO.a1}\'${DIAG.value.a1}\',${N.i5},${KD.i5},  ... ), type ${IMAT.i2}, test(${6.i2})=${RESULT[6].g12_5}';
              test.expect(RESULT[6], lessThan(THRESH), reason: reason);
              if (RESULT[6] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
              NRUN++;
            }
          }
        }, skip: skip);
      }

      // Use pathological test matrices to test DLATBS.
      for (var IMAT = NTYPE1 + 1; IMAT <= NIMAT2; IMAT++) {
        // Do the tests only if DOTYPE( IMAT ) is true.
        final skip = !DOTYPE[IMAT];

        test('DCHKTB - 2 (IN=$IN NK=$NK IMAGE=$IMAT)', () {
          final INFO = Box(0);

          for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
            // Do first for UPLO = 'U', then for UPLO = 'L'
            final UPLO = UPLOS[IUPLO - 1];
            for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
              // Do for op(A) = A, A**T, and A**H.

              final TRANS = TRANSS[ITRAN - 1];

              // Call DLATTB to generate a triangular test matrix.
              final DIAG = Box('');
              srnamc.SRNAMT = 'DLATTB';
              dlattb(IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB.asMatrix(), LDAB,
                  X, WORK, INFO);

              // +    TEST 7
              // Solve the system op(A)*x = b

              srnamc.SRNAMT = 'DLATBS';
              dcopy(N, X, 1, B, 1);
              final SCALE = Box(ZERO);
              dlatbs(UPLO, TRANS, DIAG.value, 'N', N, KD, AB.asMatrix(), LDAB,
                  B, SCALE, RWORK, INFO);

              // Check error code from DLATBS.

              if (INFO.value != 0) {
                alaerh(
                    PATH,
                    'DLATBS',
                    INFO.value,
                    0,
                    '$UPLO$TRANS${DIAG.value}N',
                    N,
                    N,
                    KD,
                    KD,
                    -1,
                    IMAT,
                    NFAIL,
                    NERRS,
                    NOUT);
              }

              dtbt03(
                  UPLO,
                  TRANS,
                  DIAG.value,
                  N,
                  KD,
                  1,
                  AB.asMatrix(),
                  LDAB,
                  SCALE.value,
                  RWORK,
                  ONE,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  WORK,
                  RESULT(7));

              // +    TEST 8
              // Solve op(A)*x = b again with NORMIN = 'Y'.

              dcopy(N, X, 1, B, 1);
              dlatbs(UPLO, TRANS, DIAG.value, 'Y', N, KD, AB.asMatrix(), LDAB,
                  B, SCALE, RWORK, INFO);

              // Check error code from DLATBS.

              if (INFO.value != 0) {
                alaerh(
                    PATH,
                    'DLATBS',
                    INFO.value,
                    0,
                    '$UPLO$TRANS${DIAG.value}Y',
                    N,
                    N,
                    KD,
                    KD,
                    -1,
                    IMAT,
                    NFAIL,
                    NERRS,
                    NOUT);
              }

              dtbt03(
                  UPLO,
                  TRANS,
                  DIAG.value,
                  N,
                  KD,
                  1,
                  AB.asMatrix(),
                  LDAB,
                  SCALE.value,
                  RWORK,
                  ONE,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  WORK,
                  RESULT(8));

              // Print information about the tests that did not pass
              // the threshold.
              for (final (K, NORMIN) in [(7, 'N'), (8, 'Y')]) {
                final reason =
                    ' DLATBS( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${NORMIN.a1}\',${N.i5},${KD.i5}, ...  ),  type ${IMAT.i2}, test(${K.i1})=${RESULT[K].g12_5}';
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(reason);
                  NFAIL++;
                }
              }
              NRUN += 2;
            }
          }
        }, skip: skip);
      }
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
