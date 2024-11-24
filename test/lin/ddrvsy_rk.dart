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
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot02.dart';
import 'dsyt01_3.dart';
import 'xlaenv.dart';

void ddrvsy_rk(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> E_,
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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final E = E_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 10, NTESTS = 3, NFACT = 2;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L']; // FACTS = ['F', 'N'];

  // Initialize constants and the random number seed.

  // Test path
  final PATH = '${'Double precision'[0]}SK';

  // Path to generate matrices
  final MATPATH = '${'Double precision'[0]}SY';

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);
  var LWORK = max(2 * NMAX, NMAX * NRHS);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrvx(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  const NB = 1;

  test.setUp(() {
    // Set the block size and minimum block size for which the block
    // routine should be used, which will be later returned by ILAENV.
    const NBMIN = 2;
    xlaenv(1, NB);
    xlaenv(2, NBMIN);
  });

  // Do for each value of N in NVAL
  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.
      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      test('DDRVSY_RK (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        // Do first for UPLO = 'U', then for UPLO = 'L'
        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];

          // Begin generate the test matrix A.

          // Set up parameters with DLATB4 for the matrix generator
          // based on the type of matrix to be generated.
          final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
              dlatb4(MATPATH, IMAT, N, N);

          // Generate a matrix with DLATMS.
          srnamc.SRNAMT = 'DLATMS';
          dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              UPLO, A.asMatrix(), LDA, WORK, INFO);

          // Check error code from DLATMS and handle error.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);

            // Skip all tests for this generated matrix
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
                var IOFF = (IZERO - 1) * LDA;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF += IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF] = ZERO;
                  IOFF += LDA;
                }
              } else {
                var IOFF = IZERO;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF] = ZERO;
                  IOFF += LDA;
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
                  IOFF += LDA;
                }
              } else {
                // Set the last IZERO rows and columns to zero.

                for (var J = 1; J <= N; J++) {
                  final I1 = max(J, IZERO);
                  for (var I = I1; I <= N; I++) {
                    A[IOFF + I] = ZERO;
                  }
                  IOFF += LDA;
                }
              }
            }
          } else {
            IZERO = 0;
          }

          // End generate the test matrix A.
          var RCONDC = ZERO;
          for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
            // Do first for FACT = 'F', then for other values.

            // Compute the condition number
            if (ZEROT) {
              if (IFACT == 1) continue;
              RCONDC = ZERO;
            } else if (IFACT == 1) {
              // Compute the 1-norm of A.
              final ANORM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);

              // Factor the matrix A.
              dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              dsytrf_rk(
                  UPLO, N, AFAC.asMatrix(), LDA, E, IWORK, WORK, LWORK, INFO);

              // Compute inv(A) and take its norm.
              dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
              LWORK = (N + NB + 1) * (NB + 3);

              // We need to compute the inverse to compute
              // RCONDC that is used later in TEST3.
              dsytri_3(
                  UPLO, N, AINV.asMatrix(), LDA, E, IWORK, WORK, LWORK, INFO);
              final AINVNM = dlansy('1', UPLO, N, AINV.asMatrix(), LDA, RWORK);

              // Compute the 1-norm condition number of A.
              if (ANORM <= ZERO || AINVNM <= ZERO) {
                RCONDC = ONE;
              } else {
                RCONDC = (ONE / ANORM) / AINVNM;
              }
            }

            // Form an exact solution and set the right hand side.
            srnamc.SRNAMT = 'DLARHS';
            dlarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';

            // --- Test DSYSV_RK  ---

            if (IFACT == 2) {
              dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              // Factor the matrix and solve the system using DSYSV_RK.
              srnamc.SRNAMT = 'DSYSV_RK';
              dsysv_rk(UPLO, N, NRHS, AFAC.asMatrix(), LDA, E, IWORK,
                  X.asMatrix(), LDA, WORK, LWORK, INFO);

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

              // Check error code from DSYSV_RK and handle error.
              test.expect(INFO.value, K);
              if (INFO.value != K) {
                alaerh(PATH, 'DSYSV_RK', INFO.value, K, UPLO, N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              } else if (INFO.value != 0) {
                //
              } else {
                // +    TEST 1      Reconstruct matrix from factors and compute
                // residual.
                dsyt01_3(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, E,
                    IWORK, AINV.asMatrix(), LDA, RWORK, RESULT(1));

                // +    TEST 2      Compute residual of the computed solution.
                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK, RESULT(2));

                // +    TEST 3
                // Check solution from generated exact solution.
                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                const NT = 3;

                // Print information about the tests that did not pass
                // the threshold.
                for (K = 1; K <= NT; K++) {
                  final reason =
                      ' DSYSV_RK, UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}';
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
          }
        }
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
