import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgesvd.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dsycon_rook.dart';
import 'package:lapack/src/dsytrf_rook.dart';
import 'package:lapack/src/dsytri_rook.dart';
import 'package:lapack/src/dsytrs_rook.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrsyx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot02.dart';
import 'dpot03.dart';
import 'dsyt01_rook.dart';
import 'xlaenv.dart';

void dchksy_rook(
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
  final Array<double> AINV_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> WORK_,
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
  const EIGHT = 8.0, SEVTEN = 17.0;
  const NTYPES = 10;
  const NTESTS = 7;
  int I, I1, I2, IOFF, IZERO, J, K, NT;
  final ISEED = Array<int>(4);
  final BLOCK = Matrix<double>(2, 2),
      DDUMMY = Array<double>(1),
      RESULT = Array<double>(NTESTS);
  final ISEEDY = Array.fromList([1988, 1989, 1990, 1991]);
  const UPLOS = ['U', 'L'];
  final INFO = Box(0), NERRS = Box(0);
  final RCOND = Box(0.0), RCONDC = Box(0.0);

  // Initialize constants and the random number seed.

  final ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  // Test path

  final PATH = '${'Double precision'[0]}SR';

  // Path to generate matrices

  final MATPATH = '${'Double precision'[0]}SY';

  var NRUN = 0;
  var NFAIL = 0;
  NERRS.value = 0;
  for (I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I];
  }

  // Test the error exits

  if (TSTERR) derrsy(PATH, NOUT);
  infoc.INFOT = 0;

  // Set the minimum block size for which the block routine should
  // be used, which will be later returned by ILAENV

  xlaenv(2, 2);

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    IZERO = 0;

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

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);

          // Skip all tests for this generated matrix

          continue;
        }

        // For matrix types 3-6, zero one or more rows and
        // columns of the matrix to test that INFO.value is returned
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
              IOFF = (IZERO - 1) * LDA;
              for (I = 1; I <= IZERO - 1; I++) {
                A[IOFF + I] = ZERO;
              }
              IOFF += IZERO;
              for (I = IZERO; I <= N; I++) {
                A[IOFF] = ZERO;
                IOFF += LDA;
              }
            } else {
              IOFF = IZERO;
              for (I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = ZERO;
                IOFF += LDA;
              }
              IOFF -= IZERO;
              for (I = IZERO; I <= N; I++) {
                A[IOFF + I] = ZERO;
              }
            }
          } else {
            if (IUPLO == 1) {
              // Set the first IZERO rows and columns to zero.

              IOFF = 0;
              for (J = 1; J <= N; J++) {
                I2 = min(J, IZERO);
                for (I = 1; I <= I2; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF += LDA;
              }
            } else {
              // Set the last IZERO rows and columns to zero.

              IOFF = 0;
              for (J = 1; J <= N; J++) {
                I1 = max(J, IZERO);
                for (I = I1; I <= N; I++) {
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
          srnamc.SRNAMT = 'DSYTRF_ROOK';
          dsytrf_rook(UPLO, N, AFAC.asMatrix(), LDA, IWORK, AINV, LWORK, INFO);

          // Adjust the expected value of INFO.value to account for
          // pivoting.

          K = IZERO;
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

          // Check error code from DSYTRF_ROOK and handle error.

          if (INFO.value != K) {
            alaerh(PATH, 'DSYTRF_ROOK', INFO.value, K, UPLO, N, N, -1, -1, NB,
                IMAT, NFAIL, NERRS, NOUT);
          }

          // Set the condition estimate flag if the INFO.value is not 0.

          final TRFCON = INFO.value != 0;

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          dsyt01_rook(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
              AINV.asMatrix(), LDA, RWORK, RESULT(1));
          NT = 1;

          // +    TEST 2
          // Form the inverse and compute the residual,
          // if the factorization was competed without INFO.value > 0
          // (i.e. there is no zero rows and columns).
          // Do it only for the first block size.

          if (INB == 1 && !TRFCON) {
            dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'DSYTRI_ROOK';
            dsytri_rook(UPLO, N, AINV.asMatrix(), LDA, IWORK, WORK, INFO);

            // Check error code from DSYTRI_ROOK and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSYTRI_ROOK', INFO.value, -1, UPLO, N, N, -1, -1,
                  -1, IMAT, NFAIL, NERRS, NOUT);
            }

            // Compute the residual for a symmetric matrix times
            // its inverse.

            dpot03(UPLO, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));
            NT = 2;
          }

          // Print information about the tests that did not pass
          // the threshold.

          for (K = 1; K <= NT; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
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

            K = N;
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

            K = 1;
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

          for (K = 3; K <= 4; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 2;

          // Skip the other tests if this is not the first block
          // size.

          if (INB > 1) continue;

          // Do only the condition estimate if INFO.value is not 0.

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

              srnamc.SRNAMT = 'DSYTRS_ROOK';
              dsytrs_rook(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK,
                  X.asMatrix(), LDA, INFO);

              // Check error code from DSYTRS_ROOK and handle error.

              if (INFO.value != 0) {
                alaerh(PATH, 'DSYTRS_ROOK', INFO.value, 0, UPLO, N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

              // Compute the residual for the solution

              dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(5));

              // +    TEST 6
              // Check solution from generated exact solution.

              dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                  RCONDC.value, RESULT(6));

              // Print information about the tests that did not pass
              // the threshold.

              for (K = 5; K <= 6; K++) {
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

          final ANORM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
          srnamc.SRNAMT = 'DSYCON_ROOK';
          dsycon_rook(UPLO, N, AFAC.asMatrix(), LDA, IWORK, ANORM, RCOND, WORK,
              IWORK(N + 1), INFO);

          // Check error code from DSYCON_ROOK and handle error.

          if (INFO.value != 0) {
            alaerh(PATH, 'DSYCON_ROOK', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                IMAT, NFAIL, NERRS, NOUT);
          }

          // Compute the test ratio to compare to values of RCOND

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
