// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

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
import 'zsyt01_rook.dart';
import 'zsyt02.dart';
import 'zsyt03.dart';

void zchksy_rook(
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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

  const ZERO = 0.0, ONE = 1.0;
  const ONEHALF = 0.5;
  const EIGHT = 8.0, SEVTEN = 17.0;
  const NTYPES = 11, NTESTS = 7;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  final BLOCK = Matrix<Complex>(2, 2), ZDUMMY = Array<Complex>(1);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  // Test path

  final PATH = '${'Zomplex precision'[0]}SR';

  // Path to generate matrices

  final MATPATH = '${'Zomplex precision'[0]}SY';

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

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Begin generate test matrix A.

        int IZERO = 0, KL = 0, KU = 0;
        if (IMAT != NTYPES) {
          // Set up parameters with ZLATB4 for the matrix generator
          // based on the type of matrix to be generated.

          final (:TYPE, KL: KLTMP, KU: KUTMP, :ANORM, :MODE, :CNDNUM, :DIST) =
              zlatb4(MATPATH, IMAT, N, N);
          (KL, KU) = (KLTMP, KUTMP);

          // Generate a matrix with ZLATMS.

          srnamc.SRNAMT = 'ZLATMS';
          zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              UPLO, A.asMatrix(), LDA, WORK, INFO);

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
                // Set the first IZERO rows and columns to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I2 = min(J, IZERO);
                  for (var I = 1; I <= I2; I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += LDA;
                }
              } else {
                // Set the last IZERO rows and columns to zero.

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

        final RCONDC = Box(ZERO);
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
          srnamc.SRNAMT = 'ZSYTRF_ROOK';
          zsytrf_rook(UPLO, N, AFAC.asMatrix(), LDA, IWORK, AINV, LWORK, INFO);

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

          // Check error code from ZSYTRF_ROOK and handle error.

          if (INFO.value != K) {
            alaerh(PATH, 'ZSYTRF_ROOK', INFO.value, K, UPLO, N, N, -1, -1, NB,
                IMAT, NFAIL, NERRS, NOUT);
          }

          // Set the condition estimate flag if the INFO is not 0.

          final TRFCON = INFO.value != 0;

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          zsyt01_rook(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
              AINV.asMatrix(), LDA, RWORK, RESULT(1));

          // +    TEST 2
          // Form the inverse and compute the residual,
          // if the factorization was competed without INFO > 0
          // (i.e. there is no zero rows and columns).
          // Do it only for the first block size.

          final int NT;
          if (INB == 1 && !TRFCON) {
            zlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'ZSYTRI_ROOK';
            zsytri_rook(UPLO, N, AINV.asMatrix(), LDA, IWORK, WORK, INFO);

            // Check error code from ZSYTRI_ROOK and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZSYTRI_ROOK', INFO.value, -1, UPLO, N, N, -1, -1,
                  -1, IMAT, NFAIL, NERRS, NOUT);
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
              NOUT.print9999(UPLO, N, NB, IMAT, K, RESULT[K]);
              NFAIL++;
            }
          }
          NRUN += NT;

          // +    TEST 3
          // Compute largest element in U or L

          RESULT[3] = ZERO;
          var DTEMP = ZERO;

          var CONST = ((pow(ALPHA, 2) - ONE) / (pow(ALPHA, 2) - ONEHALF)) /
              (ONE - ALPHA);

          if (IUPLO == 1) {
            // Compute largest element in U

            K = N;
            while (K > 1) {
              if (IWORK[K] > ZERO) {
                // Get max absolute value from elements
                // in column k in in U

                DTEMP = zlange('M', K - 1, 1,
                    AFAC((K - 1) * LDA + 1).asMatrix(), LDA, RWORK);
              } else {
                // Get max absolute value from elements
                // in columns k and k-1 in U

                DTEMP = zlange('M', K - 2, 2,
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

            K = 1;
            while (K < N) {
              if (IWORK[K] > ZERO) {
                // Get max absolute value from elements
                // in column k in in L

                DTEMP = zlange('M', N - K, 1,
                    AFAC((K - 1) * LDA + K + 1).asMatrix(), LDA, RWORK);
              } else {
                // Get max absolute value from elements
                // in columns k and k+1 in L

                DTEMP = zlange('M', N - K - 1, 2,
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
          DTEMP = ZERO;

          CONST = ((pow(ALPHA, 2) - ONE) / (pow(ALPHA, 2) - ONEHALF)) *
              ((ONE + ALPHA) / (ONE - ALPHA));

          if (IUPLO == 1) {
            // Loop backward for UPLO = 'U'

            K = N;
            while (K > 1) {
              if (IWORK[K] < ZERO) {
                // Get the two singular values
                // (real and non-negative) of a 2-by-2 block,
                // store them in RWORK array

                BLOCK[1][1] = AFAC[(K - 2) * LDA + K - 1];
                BLOCK[1][2] = AFAC[(K - 1) * LDA + K - 1];
                BLOCK[2][1] = BLOCK[1][2];
                BLOCK[2][2] = AFAC[(K - 1) * LDA + K];

                zgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, ZDUMMY.asMatrix(), 1,
                    ZDUMMY.asMatrix(), 1, WORK, 6, RWORK(3), INFO);

                final SING_MAX = RWORK[1];
                final SING_MIN = RWORK[2];

                DTEMP = SING_MAX / SING_MIN;

                // DTEMP should be bounded by CONST

                DTEMP -= CONST - THRESH;
                if (DTEMP > RESULT[4]) RESULT[4] = DTEMP;
                K--;
              }

              K--;
            }
          } else {
            // Loop forward for UPLO = 'L'

            K = 1;
            while (K < N) {
              if (IWORK[K] < ZERO) {
                // Get the two singular values
                // (real and non-negative) of a 2-by-2 block,
                // store them in RWORK array

                BLOCK[1][1] = AFAC[(K - 1) * LDA + K];
                BLOCK[2][1] = AFAC[(K - 1) * LDA + K + 1];
                BLOCK[1][2] = BLOCK[2][1];
                BLOCK[2][2] = AFAC[K * LDA + K + 1];

                zgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, ZDUMMY.asMatrix(), 1,
                    ZDUMMY.asMatrix(), 1, WORK, 6, RWORK(3), INFO);

                final SING_MAX = RWORK[1];
                final SING_MIN = RWORK[2];

                DTEMP = SING_MAX / SING_MIN;

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
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.print9999(UPLO, N, NB, IMAT, K, RESULT[K]);
              NFAIL++;
            }
          }
          NRUN += 2;

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

              // +    TEST 5 ( Using TRS_ROOK)
              // Solve and compute residual for  A * X = B.

              // Choose a set of NRHS random solution vectors
              // stored in XACT and set up the right hand side B

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(
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
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZSYTRS_ROOK';
              zsytrs_rook(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK,
                  X.asMatrix(), LDA, INFO);

              // Check error code from ZSYTRS_ROOK and handle error.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZSYTRS_ROOK', INFO.value, 0, UPLO, N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

              // Compute the residual for the solution

              zsyt02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(5));

              // +    TEST 6
              // Check solution from generated exact solution.

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                  RCONDC.value, RESULT(6));

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = 5; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(
                      ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += 2;

              // End do for each value of NRHS in NSVAL.
            }
          }

          // +    TEST 7
          // Get an estimate of RCOND = 1/CNDNUM.

          final ANORM = zlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
          final RCOND = Box(ZERO);
          srnamc.SRNAMT = 'ZSYCON_ROOK';
          zsycon_rook(
              UPLO, N, AFAC.asMatrix(), LDA, IWORK, ANORM, RCOND, WORK, INFO);

          // Check error code from ZSYCON_ROOK and handle error.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZSYCON_ROOK', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                IMAT, NFAIL, NERRS, NOUT);
          }

          // Compute the test ratio to compare values of RCOND

          RESULT[7] = dget06(RCOND.value, RCONDC.value);

          // Print information about the tests that did not pass
          // the threshold.

          if (RESULT[7] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' UPLO = \'${UPLO.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${7.i2}) =${RESULT[7].g12_5}');
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

extension on Nout {
  void print9999(String uplo, int n, int nb, int type, int test, double ratio) {
    println(
        ' UPLO = \'${uplo.a1}\', N =${n.i5}, NB =${nb.i4}, type ${type.i2}, test ${test.i2}, ratio =${ratio.g12_5}');
  }
}
