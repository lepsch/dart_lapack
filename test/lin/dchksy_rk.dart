// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrsy.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot02.dart';
import 'dpot03.dart';
import 'dsyt01_3.dart';
import 'xlaenv.dart';

void dchksy_rk(
  final Array<bool> DOTYPE_,
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
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final NSVAL = NSVAL_.having();
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
  const ZERO = 0.0, ONE = 1.0;
  const EIGHT = 8.0, SEVTEN = 17.0;
  const NTYPES = 10, NTESTS = 7;
  final BLOCK = Matrix<double>(2, 2),
      DDUMMY = Array<double>(1),
      RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];

  // Initialize constants and the random number seed.

  final ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  // Test path
  final PATH = '${'Double precision'[0]}SK';

  // Path to generate matrices
  final MATPATH = '${'Double precision'[0]}SY';

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrsy(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  test.setUp(() {
    // Set the minimum block size for which the block routine should
    // be used, which will be later returned by ILAENV
    xlaenv(2, 2);
  });

  // Do for each value of N in NVAL
  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    // Do for each value of matrix type IMAT
    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.
      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      test('DCHKSY_RK (IN=$IN IMAT=$IMAT)', () {
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

          // For matrix types 3-6, zero one or more rows and
          // columns of the matrix to test that INFO is returned
          // correctly.
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
              if (IUPLO == 1) {
                // Set the first IZERO rows and columns to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I2 = min(J, IZERO);
                  for (var I = 1; I <= I2; I++) {
                    A[IOFF + I] = ZERO;
                  }
                  IOFF += LDA;
                }
              } else {
                // Set the last IZERO rows and columns to zero.

                var IOFF = 0;
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

          // Do for each value of NB in NBVAL

          final RCONDC = Box(ZERO);
          for (var INB = 1; INB <= NNB; INB++) {
            // Set the optimal blocksize, which will be later
            // returned by ILAENV.
            final NB = NBVAL[INB];
            xlaenv(1, NB);

            // Copy the test matrix A into matrix AFAC which
            // will be factorized in place. This is needed to
            // preserve the test matrix A for subsequent tests.
            dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

            // Compute the L*D*L**T or U*D*U**T factorization of the
            // matrix. IWORK stores details of the interchanges and
            // the block structure of D. AINV is a work array for
            // block factorization, LWORK is the length of AINV.
            final LWORK = max(2, NB) * LDA;
            srnamc.SRNAMT = 'DSYTRF_RK';
            dsytrf_rk(
                UPLO, N, AFAC.asMatrix(), LDA, E, IWORK, AINV, LWORK, INFO);

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

            // Check error code from DSYTRF_RK and handle error.
            test.expect(INFO.value, K);
            if (INFO.value != K) {
              alaerh(PATH, 'DSYTRF_RK', INFO.value, K, UPLO, N, N, -1, -1, NB,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            // Set the condition estimate flag if the INFO is not 0.
            final TRFCON = INFO.value != 0;

            // +    TEST 1
            // Reconstruct matrix from factors and compute residual.

            dsyt01_3(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, E, IWORK,
                AINV.asMatrix(), LDA, RWORK, RESULT(1));

            // +    TEST 2
            // Form the inverse and compute the residual,
            // if the factorization was competed without INFO > 0
            // (i.e. there is no zero rows and columns).
            // Do it only for the first block size.

            final int NT;
            if (INB == 1 && !TRFCON) {
              dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);

              // Another reason that we need to compute the inverse
              // is that DPOT03 produces RCONDC which is used later
              // in TEST6 and TEST7.
              srnamc.SRNAMT = 'DSYTRI_3';
              final LWORK = (N + NB + 1) * (NB + 3);
              dsytri_3(
                  UPLO, N, AINV.asMatrix(), LDA, E, IWORK, WORK, LWORK, INFO);

              // Check error code from DSYTRI_3 and handle error.
              test.expect(INFO.value, 0);
              if (INFO.value != 0) {
                alaerh(PATH, 'DSYTRI_3', INFO.value, -1, UPLO, N, N, -1, -1, -1,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              // Compute the residual for a symmetric matrix times
              // its inverse.
              dpot03(UPLO, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));
              NT = 2;
            } else {
              NT = 1;
            }

            // Print information about the tests that did not pass
            // the threshold.
            for (var K = 1; K <= NT; K++) {
              final reason =
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}';
              test.expect(RESULT[K], lessThan(THRESH), reason: reason);
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += NT;

            // +    TEST 3
            // Compute largest element in U or L

            RESULT[3] = ZERO;

            var CONST = ONE / (ONE - ALPHA);

            if (IUPLO == 1) {
              // Compute largest element in U

              var K = N;
              while (K > 1) {
                double DTEMP;
                if (IWORK[K] > ZERO) {
                  // Get max absolute value from elements
                  // in column k in in U
                  DTEMP = dlange('M', K - 1, 1,
                      AFAC((K - 1) * LDA + 1).asMatrix(), LDA, RWORK);
                } else {
                  // Get max absolute value from elements
                  // in columns k and k-1 in U
                  DTEMP = dlange('M', K - 2, 2,
                      AFAC((K - 2) * LDA + 1).asMatrix(), LDA, RWORK);
                  K--;
                }

                // DTEMP should be bounded by CONST
                DTEMP -= CONST - THRESH;
                if (DTEMP > RESULT[3]) RESULT[3] = DTEMP;

                K--;
              }
            } else {
              // Compute largest element in L

              var K = 1;
              while (K < N) {
                double DTEMP;
                if (IWORK[K] > ZERO) {
                  // Get max absolute value from elements
                  // in column k in in L
                  DTEMP = dlange('M', N - K, 1,
                      AFAC((K - 1) * LDA + K + 1).asMatrix(), LDA, RWORK);
                } else {
                  // Get max absolute value from elements
                  // in columns k and k+1 in L
                  DTEMP = dlange('M', N - K - 1, 2,
                      AFAC((K - 1) * LDA + K + 2).asMatrix(), LDA, RWORK);
                  K++;
                }

                // DTEMP should be bounded by CONST
                DTEMP -= CONST - THRESH;
                if (DTEMP > RESULT[3]) RESULT[3] = DTEMP;

                K++;
              }
            }

            // +    TEST 4
            // Compute largest 2-Norm (condition number)
            // of 2-by-2 diag blocks

            RESULT[4] = ZERO;

            CONST = (ONE + ALPHA) / (ONE - ALPHA);
            dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);

            if (IUPLO == 1) {
              // Loop backward for UPLO = 'U'

              var K = N;
              while (K > 1) {
                if (IWORK[K] < ZERO) {
                  // Get the two singular values
                  // (real and non-negative) of a 2-by-2 block,
                  // store them in RWORK array

                  BLOCK[1][1] = AFAC[(K - 2) * LDA + K - 1];
                  BLOCK[1][2] = E[K];
                  BLOCK[2][1] = BLOCK[1][2];
                  BLOCK[2][2] = AFAC[(K - 1) * LDA + K];

                  dgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, DDUMMY.asMatrix(), 1,
                      DDUMMY.asMatrix(), 1, WORK, 10, INFO);

                  final SING_MAX = RWORK[1];
                  final SING_MIN = RWORK[2];

                  var DTEMP = SING_MAX / SING_MIN;

                  // DTEMP should be bounded by CONST
                  DTEMP -= CONST - THRESH;
                  if (DTEMP > RESULT[4]) RESULT[4] = DTEMP;
                  K--;
                }

                K--;
              }
            } else {
              // Loop forward for UPLO = 'L'

              var K = 1;
              while (K < N) {
                if (IWORK[K] < ZERO) {
                  // Get the two singular values
                  // (real and non-negative) of a 2-by-2 block,
                  // store them in RWORK array

                  BLOCK[1][1] = AFAC[(K - 1) * LDA + K];
                  BLOCK[2][1] = E[K];
                  BLOCK[1][2] = BLOCK[2][1];
                  BLOCK[2][2] = AFAC[K * LDA + K + 1];

                  dgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, DDUMMY.asMatrix(), 1,
                      DDUMMY.asMatrix(), 1, WORK, 10, INFO);

                  final SING_MAX = RWORK[1];
                  final SING_MIN = RWORK[2];

                  var DTEMP = SING_MAX / SING_MIN;

                  // DTEMP should be bounded by CONST

                  DTEMP -= CONST - THRESH;
                  if (DTEMP > RESULT[4]) RESULT[4] = DTEMP;
                  K++;
                }

                K++;
              }
            }

            // Print information about the tests that did not pass
            // the threshold.
            for (var K = 3; K <= 4; K++) {
              final reason =
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}';
              test.expect(RESULT[K], lessThan(THRESH), reason: reason);
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += 2;

            // Skip the other tests if this is not the first block size.
            if (INB > 1) continue;

            // Do only the condition estimate if INFO is not 0.
            if (TRFCON) {
              RCONDC.value = ZERO;
            } else {
              // Do for each value of NRHS in NSVAL.
              for (var IRHS = 1; IRHS <= NNS; IRHS++) {
                final NRHS = NSVAL[IRHS];

                // +    TEST 5 ( Using TRS_3)
                // Solve and compute residual for  A * X = B.

                // Choose a set of NRHS random solution vectors
                // stored in XACT and set up the right hand side B
                srnamc.SRNAMT = 'DLARHS';
                dlarhs(
                    MATPATH,
                    XTYPE,
                    UPLO,
                    ' ',
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
                dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                srnamc.SRNAMT = 'DSYTRS_3';
                dsytrs_3(UPLO, N, NRHS, AFAC.asMatrix(), LDA, E, IWORK,
                    X.asMatrix(), LDA, INFO);

                // Check error code from DSYTRS_3 and handle error.
                test.expect(INFO.value, 0);
                if (INFO.value != 0) {
                  alaerh(PATH, 'DSYTRS_3', INFO.value, 0, UPLO, N, N, -1, -1,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                }

                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

                // Compute the residual for the solution
                dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK, RESULT(5));

                // +    TEST 6
                // Check solution from generated exact solution.

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                    RCONDC.value, RESULT(6));

                // Print information about the tests that did not pass
                // the threshold.
                for (var K = 5; K <= 6; K++) {
                  final reason =
                      ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}';
                  test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.println();
                    NFAIL++;
                  }
                }
                NRUN += 2;
                // End do for each value of NRHS in NSVAL.
              }
            }

            // +    TEST 7
            // Get an estimate of RCOND = 1/CNDNUM.

            final ANORM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
            final RCOND = Box(0.0);
            srnamc.SRNAMT = 'DSYCON_3';
            dsycon_3(UPLO, N, AFAC.asMatrix(), LDA, E, IWORK, ANORM, RCOND,
                WORK, IWORK(N + 1), INFO);

            // Check error code from DSYCON_3 and handle error.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DSYCON_3', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            // Compute the test ratio to compare to values of RCOND
            RESULT[7] = dget06(RCOND.value, RCONDC.value);

            // Print information about the tests that did not pass
            // the threshold.
            final reason =
                ' UPLO = \'${UPLO.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${7.i2}) =${RESULT[7].g12_5}';
            test.expect(RESULT[7], lessThan(THRESH), reason: reason);
            if (RESULT[7] >= THRESH) {
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

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
