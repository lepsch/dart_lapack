// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrge.dart';
import 'dget01.dart';
import 'dget02.dart';
import 'dget03.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dget07.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'xlaenv.dart';

void dchkge(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
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
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 11, NTESTS = 8, NTRAN = 3;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991], TRANSS = ['N', 'T', 'C'];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}GE';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  // Test the error exits

  test.group('error exits', () {
    test.setUp(() {
      xlaenv(1, 1);
    });
    if (TSTERR) derrge(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  test.setUp(() {
    xlaenv(2, 2);
  });

  // Do for each value of M in MVAL
  for (final IM in 1.through(NM)) {
    final M = MVAL[IM];
    final LDA = max(1, M);

    // Do for each value of N in NVAL
    for (final IN in 1.through(NN)) {
      final N = NVAL[IN];
      final NIMAT = M <= 0 || N <= 0 ? 1 : NTYPES;

      for (final IMAT in 1.through(NIMAT)) {
        // Do the tests only if DOTYPE( IMAT ) is true.
        final skip = !DOTYPE[IMAT];

        // Skip types 5, 6, or 7 if the matrix size is too small.
        final ZEROT = IMAT >= 5 && IMAT <= 7;
        if (ZEROT && N < IMAT - 4) continue;

        test('DCHKGE (IM=$IM IN=$IN IMAT=$IMAT)', () {
          final INFO = Box(0);

          // Set up parameters with DLATB4 and generate a test matrix
          // with DLATMS.

          final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
              dlatb4(PATH, IMAT, M, N);

          srnamc.SRNAMT = 'DLATMS';
          dlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              'No packing', A.asMatrix(), LDA, WORK, INFO);

          // Check error code from DLATMS.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            return;
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
            final IOFF = (IZERO - 1) * LDA;
            if (IMAT < 7) {
              for (var I = 1; I <= M; I++) {
                A[IOFF + I] = ZERO;
              }
            } else {
              dlaset('Full', M, N - IZERO + 1, ZERO, ZERO,
                  A(IOFF + 1).asMatrix(), LDA);
            }
          } else {
            IZERO = 0;
          }

          // These lines, if used in place of the calls in the DO 60
          // loop, cause the code to bomb on a Sun SPARCstation.

          // ANORMO = dlange( 'O', M, N, A, LDA, RWORK )
          // ANORMI = dlange( 'I', M, N, A, LDA, RWORK )

          // Do for each blocksize in NBVAL

          for (var INB = 1; INB <= NNB; INB++) {
            final NB = NBVAL[INB];
            xlaenv(1, NB);

            // Compute the LU factorization of the matrix.

            dlacpy('Full', M, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            srnamc.SRNAMT = 'DGETRF';
            dgetrf(M, N, AFAC.asMatrix(), LDA, IWORK, INFO);

            // Check error code from DGETRF.
            test.expect(INFO.value, IZERO);
            if (INFO.value != IZERO) {
              alaerh(PATH, 'DGETRF', INFO.value, IZERO, ' ', M, N, -1, -1, NB,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            // TEST 1
            // Reconstruct matrix from factors and compute residual.

            dlacpy('Full', M, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            dget01(M, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA, IWORK, RWORK,
                RESULT(1));

            // TEST 2
            // Form the inverse if the factorization was successful
            // and compute the residual.

            final int NT;
            final bool TRFCON;
            final double ANORMI, ANORMO, RCONDI;
            final RCONDO = Box(0.0);
            if (M == N && INFO.value == 0) {
              TRFCON = false;
              dlacpy('Full', N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
              srnamc.SRNAMT = 'DGETRI';
              final NRHS = NSVAL[1];
              final LWORK = NMAX * max(3, NRHS).toInt();
              dgetri(N, AINV.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

              // Check error code from DGETRI.
              test.expect(INFO.value, 0, reason: 'DGETRI');
              if (INFO.value != 0) {
                alaerh(PATH, 'DGETRI', INFO.value, 0, ' ', N, N, -1, -1, NB,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              // Compute the residual for the matrix times its
              // inverse.  Also compute the 1-norm condition number
              // of A.

              dget03(N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RCONDO, RESULT(2));
              ANORMO = dlange('O', M, N, A.asMatrix(), LDA, RWORK);

              // Compute the infinity-norm condition number of A.

              ANORMI = dlange('I', M, N, A.asMatrix(), LDA, RWORK);
              final AINVNM = dlange('I', N, N, AINV.asMatrix(), LDA, RWORK);
              if (ANORMI <= ZERO || AINVNM <= ZERO) {
                RCONDI = ONE;
              } else {
                RCONDI = (ONE / ANORMI) / AINVNM;
              }
              NT = 2;
            } else {
              // Do only the condition estimate if INFO > 0.

              TRFCON = true;
              ANORMO = dlange('O', M, N, A.asMatrix(), LDA, RWORK);
              ANORMI = dlange('I', M, N, A.asMatrix(), LDA, RWORK);
              RCONDO.value = ZERO;
              RCONDI = ZERO;
              NT = 1;
            }

            // Print information about the tests so far that did not
            // pass the threshold.

            for (var K = 1; K <= NT; K++) {
              final reason =
                  ' M = ${M.i5}, N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}';
              test.expect(RESULT[K], lessThan(THRESH), reason: reason);
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += NT;

            // Skip the remaining tests if this is not the first
            // block size or if M != N.  Skip the solve tests if
            // the matrix is singular.

            if (INB > 1 || M != N) continue;
            if (!TRFCON) {
              for (var IRHS = 1; IRHS <= NNS; IRHS++) {
                final NRHS = NSVAL[IRHS];
                var XTYPE = 'N';

                for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
                  final TRANS = TRANSS[ITRAN - 1];
                  final RCONDC = ITRAN == 1 ? RCONDO.value : RCONDI;

                  // TEST 3
                  // Solve and compute residual for A * X = B.

                  srnamc.SRNAMT = 'DLARHS';
                  dlarhs(
                      PATH,
                      XTYPE,
                      ' ',
                      TRANS,
                      N,
                      N,
                      KL,
                      KU,
                      NRHS,
                      A.asMatrix(),
                      LDA,
                      XACT.asMatrix(),
                      LDA,
                      B.asMatrix(),
                      LDA,
                      ISEED,
                      INFO);
                  XTYPE = 'C';

                  dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);
                  srnamc.SRNAMT = 'DGETRS';
                  dgetrs(TRANS, N, NRHS, AFAC.asMatrix(), LDA, IWORK,
                      X.asMatrix(), LDA, INFO);

                  // Check error code from DGETRS.
                  test.expect(INFO.value, 0, reason: 'DGETRS');
                  if (INFO.value != 0) {
                    alaerh(PATH, 'DGETRS', INFO.value, 0, TRANS, N, N, -1, -1,
                        NRHS, IMAT, NFAIL, NERRS, NOUT);
                  }

                  dlacpy(
                      'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                  dget02(TRANS, N, N, NRHS, A.asMatrix(), LDA, X.asMatrix(),
                      LDA, WORK.asMatrix(), LDA, RWORK, RESULT(3));

                  // TEST 4
                  // Check solution from generated exact solution.

                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(4));

                  // TESTS 5, 6, and 7
                  // Use iterative refinement to improve the
                  // solution.

                  srnamc.SRNAMT = 'DGERFS';
                  dgerfs(
                      TRANS,
                      N,
                      NRHS,
                      A.asMatrix(),
                      LDA,
                      AFAC.asMatrix(),
                      LDA,
                      IWORK,
                      B.asMatrix(),
                      LDA,
                      X.asMatrix(),
                      LDA,
                      RWORK,
                      RWORK(NRHS + 1),
                      WORK,
                      IWORK(N + 1),
                      INFO);

                  // Check error code from DGERFS.
                  test.expect(INFO.value, 0, reason: 'DGERFS');
                  if (INFO.value != 0) {
                    alaerh(PATH, 'DGERFS', INFO.value, 0, TRANS, N, N, -1, -1,
                        NRHS, IMAT, NFAIL, NERRS, NOUT);
                  }

                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(5));
                  dget07(
                      TRANS,
                      N,
                      NRHS,
                      A.asMatrix(),
                      LDA,
                      B.asMatrix(),
                      LDA,
                      X.asMatrix(),
                      LDA,
                      XACT.asMatrix(),
                      LDA,
                      RWORK,
                      true,
                      RWORK(NRHS + 1),
                      RESULT(6));

                  // Print information about the tests that did not
                  // pass the threshold.

                  for (var K = 3; K <= 7; K++) {
                    final reason =
                        ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}';
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
            }

            // TEST 8
            // Get an estimate of RCOND = 1/CNDNUM.

            for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
              final String NORM;
              final double ANORM, RCONDC;
              if (ITRAN == 1) {
                ANORM = ANORMO;
                RCONDC = RCONDO.value;
                NORM = 'O';
              } else {
                ANORM = ANORMI;
                RCONDC = RCONDI;
                NORM = 'I';
              }
              srnamc.SRNAMT = 'DGECON';
              final RCOND = Box(0.0);
              dgecon(NORM, N, AFAC.asMatrix(), LDA, ANORM, RCOND, WORK,
                  IWORK(N + 1), INFO);

              // Check error code from DGECON.
              test.expect(INFO.value, 0, reason: 'DGECON');
              if (INFO.value != 0) {
                alaerh(PATH, 'DGECON', INFO.value, 0, NORM, N, N, -1, -1, -1,
                    IMAT, NFAIL, NERRS, NOUT);
              }
              RESULT[8] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.
              final reason =
                  ' NORM =\'${NORM.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${8.i2}) =${RESULT[8].g12_5}';
              test.expect(RESULT[8], lessThan(THRESH), reason: reason);
              if (RESULT[8] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
              NRUN++;
            }
          }
        }, skip: skip);
      }
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
