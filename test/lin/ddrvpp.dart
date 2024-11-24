// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';
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
import 'dppt01.dart';
import 'dppt02.dart';
import 'dppt05.dart';

void ddrvpp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> ASAV_,
  final Array<double> B_,
  final Array<double> BSAV_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> S_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final ASAV = ASAV_.having();
  final B = B_.having();
  final BSAV = BSAV_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 9, NTESTS = 6;
  double ROLDC = 0;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'],
      FACTS = ['F', 'N', 'E'],
      PACKS = ['C', 'R'],
      EQUEDS = ['N', 'Y'];
  final AMAX = Box(ZERO), SCOND = Box(ZERO), RCOND = Box(ZERO);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}PP';
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
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      // Skip types 3, 4, or 5 if the matrix size is too small.
      final ZEROT = IMAT >= 3 && IMAT <= 5;
      if (ZEROT && N < IMAT - 2) continue;

      test('DDRVPP (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        // Do first for UPLO = 'U', then for UPLO = 'L'
        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];
          final PACKIT = PACKS[IUPLO - 1];

          // Set up parameters with DLATB4 and generate a test matrix
          // with DLATMS.

          final (:TYPE, :KL, :KU, :ANORM, :MODE, :COND, :DIST) =
              dlatb4(PATH, IMAT, N, N);
          var RCONDC = ONE / COND;

          srnamc.SRNAMT = 'DLATMS';
          dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU,
              PACKIT, A.asMatrix(), LDA, WORK, INFO);

          // Check error code from DLATMS.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            continue;
          }

          // For types 3-5, zero one row and column of the matrix to
          // test that INFO is returned correctly.
          final int IZERO;
          if (ZEROT) {
            if (IMAT == 3) {
              IZERO = 1;
            } else if (IMAT == 4) {
              IZERO = N;
            } else {
              IZERO = N ~/ 2 + 1;
            }

            // Set row and column IZERO of A to 0.

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
            IZERO = 0;
          }

          // Save a copy of the matrix A in ASAV.
          dcopy(NPP, A, 1, ASAV, 1);

          for (var IEQUED = 1; IEQUED <= 2; IEQUED++) {
            final EQUED = Box(EQUEDS[IEQUED - 1]);
            final NFACT = IEQUED == 1 ? 3 : 1;

            for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
              final FACT = FACTS[IFACT - 1];
              final PREFAC = lsame(FACT, 'F');
              final NOFACT = lsame(FACT, 'N');
              final EQUIL = lsame(FACT, 'E');

              if (ZEROT) {
                if (PREFAC) continue;
                RCONDC = ZERO;
              } else if (!lsame(FACT, 'N')) {
                // Compute the condition number for comparison with
                // the value returned by DPPSVX (FACT = 'N' reuses
                // the condition number from the previous iteration
                // with FACT = 'F').

                dcopy(NPP, ASAV, 1, AFAC, 1);
                if (EQUIL || IEQUED > 1) {
                  // Compute row and column scale factors to
                  // equilibrate the matrix A.

                  dppequ(UPLO, N, AFAC, S, SCOND, AMAX, INFO);
                  if (INFO.value == 0 && N > 0) {
                    if (IEQUED > 1) SCOND.value = ZERO;

                    // Equilibrate the matrix.

                    dlaqsp(UPLO, N, AFAC, S, SCOND.value, AMAX.value, EQUED);
                  }
                }

                // Save the condition number of the
                // non-equilibrated system for use in DGET04.

                if (EQUIL) ROLDC = RCONDC;

                // Compute the 1-norm of A.

                final ANORM = dlansp('1', UPLO, N, AFAC, RWORK);

                // Factor the matrix A.

                dpptrf(UPLO, N, AFAC, INFO);

                // Form the inverse of A.

                dcopy(NPP, AFAC, 1, A, 1);
                dpptri(UPLO, N, A, INFO);

                // Compute the 1-norm condition number of A.

                final AINVNM = dlansp('1', UPLO, N, A, RWORK);
                if (ANORM <= ZERO || AINVNM <= ZERO) {
                  RCONDC = ONE;
                } else {
                  RCONDC = (ONE / ANORM) / AINVNM;
                }
              }

              // Restore the matrix A.

              dcopy(NPP, ASAV, 1, A, 1);

              // Form an exact solution and set the right hand side.

              srnamc.SRNAMT = 'DLARHS';
              dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                  LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
              XTYPE = 'C';
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

              if (NOFACT) {
                // --- Test DPPSV  ---

                // Compute the L*L' or U'*U factorization of the
                // matrix and solve the system.

                dcopy(NPP, A, 1, AFAC, 1);
                dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                srnamc.SRNAMT = 'DPPSV';
                dppsv(UPLO, N, NRHS, AFAC, X.asMatrix(), LDA, INFO);

                // Check error code from DPPSV .
                test.expect(INFO.value, IZERO);
                if (INFO.value != IZERO) {
                  alaerh(PATH, 'DPPSV ', INFO.value, IZERO, UPLO, N, N, -1, -1,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                  continue;
                } else if (INFO.value != 0) {
                  continue;
                }

                // Reconstruct matrix from factors and compute
                // residual.

                dppt01(UPLO, N, A, AFAC, RWORK, RESULT(1));

                // Compute residual of the computed solution.

                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(),
                    LDA, RWORK, RESULT(2));

                // Check solution from generated exact solution.

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                const NT = 3;

                // Print information about the tests that did not
                // pass the threshold.

                for (var K = 1; K <= NT; K++) {
                  final reason =
                      ' DPPSV, UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}';
                  test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                    NOUT.println(reason);
                    NFAIL++;
                  }
                }
                NRUN += NT;
              }

              // --- Test DPPSVX ---

              if (!PREFAC && NPP > 0) {
                dlaset('Full', NPP, 1, ZERO, ZERO, AFAC.asMatrix(), NPP);
              }
              dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT='F' and
                // EQUED='Y'.

                dlaqsp(UPLO, N, A, S, SCOND.value, AMAX.value, EQUED);
              }

              // Solve the system and compute the condition number
              // and error bounds using DPPSVX.

              srnamc.SRNAMT = 'DPPSVX';
              dppsvx(
                  FACT,
                  UPLO,
                  N,
                  NRHS,
                  A,
                  AFAC,
                  EQUED,
                  S,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  RCOND,
                  RWORK,
                  RWORK(NRHS + 1),
                  WORK,
                  IWORK,
                  INFO);

              // Check the error code from DPPSVX.
              test.expect(INFO.value, IZERO);
              if (INFO.value != IZERO) {
                alaerh(PATH, 'DPPSVX', INFO.value, IZERO, FACT + UPLO, N, N, -1,
                    -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }
              final int K1;
              if (INFO.value == 0) {
                if (!PREFAC) {
                  // Reconstruct matrix from factors and compute
                  // residual.

                  dppt01(UPLO, N, A, AFAC, RWORK(2 * NRHS + 1), RESULT(1));
                  K1 = 1;
                } else {
                  K1 = 2;
                }

                // Compute residual of the computed solution.

                dlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                dppt02(UPLO, N, NRHS, ASAV, X.asMatrix(), LDA, WORK.asMatrix(),
                    LDA, RWORK(2 * NRHS + 1), RESULT(2));

                // Check solution from generated exact solution.

                if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                } else {
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      ROLDC, RESULT(3));
                }

                // Check the error bounds from iterative
                // refinement.

                dppt05(
                    UPLO,
                    N,
                    NRHS,
                    ASAV,
                    B.asMatrix(),
                    LDA,
                    X.asMatrix(),
                    LDA,
                    XACT.asMatrix(),
                    LDA,
                    RWORK,
                    RWORK(NRHS + 1),
                    RESULT(4));
              } else {
                K1 = 6;
              }

              // Compare RCOND from DPPSVX with the computed value
              // in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = K1; K <= 6; K++) {
                final reason = PREFAC
                    ? ' DPPSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}'
                    : ' DPPSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, EQUED=\'${EQUED.value.a1}\', type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}';
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
        }
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
