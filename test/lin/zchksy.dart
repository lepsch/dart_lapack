// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zerrsy.dart';
import 'zget04.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zlatsy.dart';
import 'zpot05.dart';
import 'zsyt01.dart';
import 'zsyt02.dart';
import 'zsyt03.dart';

void zchksy(
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
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<Complex> AINV_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
  final DOTYPE = DOTYPE_.having();
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
  const ZERO = 0.0;
  const NTYPES = 11, NTESTS = 9;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}SY';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrsy(PATH, NOUT);
  infoc.INFOT = 0;

  // Set the minimum block size for which the block routine should
  // be used, which will be later returned by ILAENV

  xlaenv(2, 2);

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    // Do for each value of matrix type IMAT

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      var IZERO = 0, KL = 0, KU = 0;
      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Begin generate test matrix A.

        if (IMAT != NTYPES) {
          // Set up parameters with ZLATB4 for the matrix generator
          // based on the type of matrix to be generated.

          final (:TYPE, KL: KLTMP, KU: KUTMP, :ANORM, :MODE, :CNDNUM, :DIST) =
              zlatb4(PATH, IMAT, N, N);
          (KL, KU) = (KLTMP, KUTMP);

          // Generate a matrix with ZLATMS.

          srnamc.SRNAMT = 'ZLATMS';
          zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              'N', A.asMatrix(), LDA, WORK, INFO);

          // Check error code from ZLATMS and handle error.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);

            // Skip all tests for this generated matrix

            continue;
          }

          // For matrix types 3-6, zero one or more rows and
          // columns of the matrix to test that INFO is returned
          // correctly.

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
                  A[IOFF + I] = Complex.zero;
                }
                IOFF += IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF] = Complex.zero;
                  IOFF += LDA;
                }
              } else {
                var IOFF = IZERO;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF] = Complex.zero;
                  IOFF += LDA;
                }
                IOFF -= IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF + I] = Complex.zero;
                }
              }
            } else {
              if (IUPLO == 1) {
                // Set the first IZERO rows to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I2 = min(J, IZERO);
                  for (var I = 1; I <= I2; I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += LDA;
                }
              } else {
                // Set the last IZERO rows to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I1 = max(J, IZERO);
                  for (var I = I1; I <= N; I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += LDA;
                }
              }
            }
          } else {
            IZERO = 0;
          }
        } else {
          // For matrix kind IMAT = 11, generate special block
          // diagonal matrix to test alternate code
          // for the 2 x 2 blocks.

          zlatsy(UPLO, N, A.asMatrix(), LDA, ISEED);
        }

        // End generate test matrix A.

        // Do for each value of NB in NBVAL

        for (var INB = 1; INB <= NNB; INB++) {
          // Set the optimal blocksize, which will be later
          // returned by ILAENV.

          final NB = NBVAL[INB];
          xlaenv(1, NB);

          // Copy the test matrix A into matrix AFAC which
          // will be factorized in place. This is needed to
          // preserve the test matrix A for subsequent tests.

          zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

          // Compute the L*D*L**T or U*D*U**T factorization of the
          // matrix. IWORK stores details of the interchanges and
          // the block structure of D. AINV is a work array for
          // block factorization, LWORK is the length of AINV.

          final LWORK = max(2, NB) * LDA;
          srnamc.SRNAMT = 'ZSYTRF';
          zsytrf(UPLO, N, AFAC.asMatrix(), LDA, IWORK, AINV, LWORK, INFO);

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

          // Check error code from ZSYTRF and handle error.

          if (INFO.value != K) {
            alaerh(PATH, 'ZSYTRF', INFO.value, K, UPLO, N, N, -1, -1, NB, IMAT,
                NFAIL, NERRS, NOUT);
          }

          // Set the condition estimate flag if the INFO is not 0.

          final TRFCON = INFO.value != 0;

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          zsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
              AINV.asMatrix(), LDA, RWORK, RESULT(1));

          // +    TEST 2
          // Form the inverse and compute the residual,
          // if the factorization was competed without INFO > 0
          // (i.e. there is no zero rows and columns).
          // Do it only for the first block size.

          final int NT;
          final RCONDC = Box(ZERO);
          if (INB == 1 && !TRFCON) {
            zlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'ZSYTRI2';
            final LWORK = (N + NB + 1) * (NB + 3);
            zsytri2(UPLO, N, AINV.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Check error code from ZSYTRI2 and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZSYTRI2', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            // Compute the residual for a symmetric matrix times
            // its inverse.

            zsyt03(UPLO, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));
            NT = 2;
          } else {
            NT = 1;
          }

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NT; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += NT;

          // Skip the other tests if this is not the first block
          // size.

          if (INB > 1) continue;

          // Do only the condition estimate if INFO is not 0.

          if (TRFCON) {
            RCONDC.value = ZERO;
          } else {
            // Do for each value of NRHS in NSVAL.

            for (var IRHS = 1; IRHS <= NNS; IRHS++) {
              final NRHS = NSVAL[IRHS];

              // +    TEST 3 (Using TRS)
              // Solve and compute residual for  A * X = B.

              // Choose a set of NRHS random solution vectors
              // stored in XACT and set up the right hand side B

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                  LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZSYTRS';
              zsytrs(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(),
                  LDA, INFO);

              // Check error code from ZSYTRS and handle error.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZSYTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

              // Compute the residual for the solution

              zsyt02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(3));

              // +    TEST 4 (Using TRS2)
              // Solve and compute residual for  A * X = B.

              // Choose a set of NRHS random solution vectors
              // stored in XACT and set up the right hand side B

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                  LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZSYTRS2';
              zsytrs2(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(),
                  LDA, WORK, INFO);

              // Check error code from ZSYTRS2 and handle error.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZSYTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

              // Compute the residual for the solution

              zsyt02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(4));

              // +    TEST 5
              // Check solution from generated exact solution.

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                  RCONDC.value, RESULT(5));

              // +    TESTS 6, 7, and 8
              // Use iterative refinement to improve the solution.

              srnamc.SRNAMT = 'ZSYRFS';
              zsyrfs(
                  UPLO,
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
                  RWORK(2 * NRHS + 1),
                  INFO);

              // Check error code from ZSYRFS and handle error.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZSYRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                  RCONDC.value, RESULT(6));
              zpot05(
                  UPLO,
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
                  RWORK(NRHS + 1),
                  RESULT(7));

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = 3; K <= 8; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(
                      ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += 6;

              // End do for each value of NRHS in NSVAL.
            }
          }

          // +    TEST 9
          // Get an estimate of RCOND = 1/CNDNUM.

          final ANORM = zlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
          srnamc.SRNAMT = 'ZSYCON';
          final RCOND = Box(ZERO);
          zsycon(
              UPLO, N, AFAC.asMatrix(), LDA, IWORK, ANORM, RCOND, WORK, INFO);

          // Check error code from ZSYCON and handle error.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZSYCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          // Compute the test ratio to compare values of RCOND

          RESULT[9] = dget06(RCOND.value, RCONDC.value);

          // Print information about the tests that did not pass
          // the threshold.

          if (RESULT[9] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' UPLO = \'${UPLO.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${9.i2}) =${RESULT[9].g12_5}');
            NFAIL++;
          }
          NRUN++;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
