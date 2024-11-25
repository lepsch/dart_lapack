// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dppt02.dart';
import 'dppt05.dart';
import 'dspt01.dart';

void ddrvsp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
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
  final NVAL = NVAL_.having();
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
  const NTYPES = 10, NTESTS = 6, NFACT = 2;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const FACTS = ['F', 'N'];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}SP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrvx(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  // Do for each value of N in NVAL
  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final NPP = N * (N + 1) ~/ 2;
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.
      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      test('DDRVSP (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);
        String? XTYPE;

        // Do first for UPLO = 'U', then for UPLO = 'L'
        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final (UPLO, PACKIT) = IUPLO == 1 ? ('U', 'C') : ('L', 'R');

          // Set up parameters with DLATB4 and generate a test matrix
          // with DLATMS.
          final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
              dlatb4(PATH, IMAT, N, N);

          srnamc.SRNAMT = 'DLATMS';
          dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              PACKIT, A.asMatrix(), LDA, WORK, INFO);

          // Check error code from DLATMS.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            continue;
          }

          // For types 3-6, zero one or more rows and columns of the
          // matrix to test that INFO is returned correctly.

          final int IZERO;
          if (ZEROT) {
            if (IMAT == 3) {
              IZERO = 1;
            } else if (IMAT == 4) {
              IZERO = N;
            } else {
              IZERO = N ~/ 2 + 1;
            }

            if (IMAT < 6) {
              // Set row and column IZERO to zero.

              if (IUPLO == 1) {
                var IOFF = (IZERO - 1) * IZERO ~/ 2;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF += IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF] = ZERO;
                  IOFF += I;
                }
              } else {
                var IOFF = IZERO;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF] = ZERO;
                  IOFF += N - I;
                }
                IOFF -= IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF + I] = ZERO;
                }
              }
            } else {
              var IOFF = 0;
              if (IUPLO == 1) {
                // Set the first IZERO rows and columns to zero.

                for (var J = 1; J <= N; J++) {
                  final I2 = min(J, IZERO);
                  for (var I = 1; I <= I2; I++) {
                    A[IOFF + I] = ZERO;
                  }
                  IOFF += J;
                }
              } else {
                // Set the last IZERO rows and columns to zero.

                for (var J = 1; J <= N; J++) {
                  final I1 = max(J, IZERO);
                  for (var I = I1; I <= N; I++) {
                    A[IOFF + I] = ZERO;
                  }
                  IOFF += N - J;
                }
              }
            }
          } else {
            IZERO = 0;
          }

          var RCONDC = 0.0;
          for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
            // Do first for FACT = 'F', then for other values.
            final FACT = FACTS[IFACT - 1];

            // Compute the condition number for comparison with
            // the value returned by DSPSVX.

            if (ZEROT) {
              if (IFACT == 1) continue;
              RCONDC = ZERO;
            } else if (IFACT == 1) {
              // Compute the 1-norm of A.
              final ANORM = dlansp('1', UPLO, N, A, RWORK);

              // Factor the matrix A.
              dcopy(NPP, A, 1, AFAC, 1);
              dsptrf(UPLO, N, AFAC, IWORK, INFO);

              // Compute inv(A) and take its norm.
              dcopy(NPP, AFAC, 1, AINV, 1);
              dsptri(UPLO, N, AINV, IWORK, WORK, INFO);
              final AINVNM = dlansp('1', UPLO, N, AINV, RWORK);

              // Compute the 1-norm condition number of A.
              if (ANORM <= ZERO || AINVNM <= ZERO) {
                RCONDC = ONE;
              } else {
                RCONDC = (ONE / ANORM) / AINVNM;
              }
            }

            // Form an exact solution and set the right hand side.
            srnamc.SRNAMT = 'DLARHS';
            XTYPE ??= IMAT == 1 ? 'N' : 'C';
            dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';

            // --- Test DSPSV  ---

            if (IFACT == 2) {
              dcopy(NPP, A, 1, AFAC, 1);
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              // Factor the matrix and solve the system using DSPSV.
              srnamc.SRNAMT = 'DSPSV';
              dspsv(UPLO, N, NRHS, AFAC, IWORK, X.asMatrix(), LDA, INFO);

              // Adjust the expected value of INFO to account for
              // pivoting.
              var K = IZERO;
              if (K > 0) {
                while (true) {
                  if (IWORK[K] < 0) {
                    if (IWORK[K] != -K) {
                      K = -IWORK[K];
                      continue;
                    }
                  } else if (IWORK[K] != K) {
                    K = IWORK[K];
                    continue;
                  }
                  break;
                }
              }

              // Check error code from DSPSV .
              test.expect(INFO.value, K);
              if (INFO.value != K) {
                alaerh(PATH, 'DSPSV ', INFO.value, K, UPLO, N, N, -1, -1, NRHS,
                    IMAT, NFAIL, NERRS, NOUT);
              } else if (INFO.value != 0) {
                //
              } else {
                // Reconstruct matrix from factors and compute residual.
                dspt01(UPLO, N, A, AFAC, IWORK, AINV.asMatrix(), LDA, RWORK,
                    RESULT(1));

                // Compute residual of the computed solution.
                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(),
                    LDA, RWORK, RESULT(2));

                // Check solution from generated exact solution.
                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                const NT = 3;

                // Print information about the tests that did not pass
                // the threshold.
                for (var K = 1; K <= NT; K++) {
                  final reason =
                      ' DSPSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}';
                  test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                    NOUT.println(reason);
                    NFAIL++;
                  }
                }
                NRUN += NT;
              }
            }

            // --- Test DSPSVX ---

            if (IFACT == 2 && NPP > 0) {
              dlaset('Full', NPP, 1, ZERO, ZERO, AFAC.asMatrix(), NPP);
            }
            dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);

            // Solve the system and compute the condition number and
            // error bounds using DSPSVX.
            final RCOND = Box(ZERO);
            srnamc.SRNAMT = 'DSPSVX';
            dspsvx(
                FACT,
                UPLO,
                N,
                NRHS,
                A,
                AFAC,
                IWORK,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                RCOND,
                RWORK,
                RWORK(NRHS + 1),
                WORK,
                IWORK(N + 1),
                INFO);

            // Adjust the expected value of INFO to account for pivoting.
            var K = IZERO;
            if (K > 0) {
              while (true) {
                if (IWORK[K] < 0) {
                  if (IWORK[K] != -K) {
                    K = -IWORK[K];
                    continue;
                  }
                } else if (IWORK[K] != K) {
                  K = IWORK[K];
                  continue;
                }
                break;
              }
            }

            // Check the error code from DSPSVX.
            test.expect(INFO.value, K);
            if (INFO.value != K) {
              alaerh(PATH, 'DSPSVX', INFO.value, K, FACT + UPLO, N, N, -1, -1,
                  NRHS, IMAT, NFAIL, NERRS, NOUT);
              continue;
            }

            final int K1;
            if (INFO.value == 0) {
              if (IFACT >= 2) {
                // Reconstruct matrix from factors and compute residual.
                dspt01(UPLO, N, A, AFAC, IWORK, AINV.asMatrix(), LDA,
                    RWORK(2 * NRHS + 1), RESULT(1));
                K1 = 1;
              } else {
                K1 = 2;
              }

              // Compute residual of the computed solution.
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                  RWORK(2 * NRHS + 1), RESULT(2));

              // Check solution from generated exact solution.
              dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));

              // Check the error bounds from iterative refinement.
              dppt05(UPLO, N, NRHS, A, B.asMatrix(), LDA, X.asMatrix(), LDA,
                  XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(4));
            } else {
              K1 = 6;
            }

            // Compare RCOND from DSPSVX with the computed value in RCONDC.
            RESULT[6] = dget06(RCOND.value, RCONDC);

            // Print information about the tests that did not pass
            // the threshold.
            for (var K = K1; K <= 6; K++) {
              final reason =
                  ' DSPSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}';
              test.expect(RESULT[K], lessThan(THRESH), reason: reason);
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += 7 - K1;
          }
        }
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
