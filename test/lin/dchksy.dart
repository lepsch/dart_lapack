import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dsycon.dart';
import 'package:lapack/src/dsyrfs.dart';
import 'package:lapack/src/dsytrf.dart';
import 'package:lapack/src/dsytri2.dart';
import 'package:lapack/src/dsytrs.dart';
import 'package:lapack/src/dsytrs2.dart';
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
import 'dpot05.dart';
import 'dsyt01.dart';
import 'xlaenv.dart';

void dchksy(
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

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const NTYPES = 10, NTESTS = 9;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);
  final RCONDC = Box(0.0), RCOND = Box(0.0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}SY';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
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
            dlatb4(PATH, IMAT, N, N);

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
              IOFF = IOFF + IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF] = ZERO;
                IOFF = IOFF + LDA;
              }
            } else {
              var IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = ZERO;
                IOFF = IOFF + LDA;
              }
              IOFF = IOFF - IZERO;
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
                IOFF = IOFF + LDA;
              }
            } else {
              // Set the last IZERO rows and columns to zero.

              var IOFF = 0;
              for (var J = 1; J <= N; J++) {
                var I1 = max(J, IZERO);
                for (var I = I1; I <= N; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF = IOFF + LDA;
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
          srnamc.SRNAMT = 'DSYTRF';
          dsytrf(UPLO, N, AFAC.asMatrix(), LDA, IWORK, AINV, LWORK, INFO);

          // Adjust the expected value of INFO.value to account for
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

          // Check error code from DSYTRF and handle error.

          if (INFO.value != K) {
            alaerh(PATH, 'DSYTRF', INFO.value, K, UPLO, N, N, -1, -1, NB, IMAT,
                NFAIL, NERRS, NOUT);
          }

          // Set the condition estimate flag if the INFO is not 0.

          final TRFCON = INFO.value != 0;

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          dsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
              AINV.asMatrix(), LDA, RWORK, RESULT(1));
          var NT = 1;

          // +    TEST 2
          // Form the inverse and compute the residual,
          // if the factorization was competed without INFO.value > 0
          // (i.e. there is no zero rows and columns).
          // Do it only for the first block size.

          if (INB == 1 && !TRFCON) {
            dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'DSYTRI2';
            final LWORK = (N + NB + 1) * (NB + 3);
            dsytri2(UPLO, N, AINV.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Check error code from DSYTRI2 and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSYTRI2', INFO.value, -1, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            // Compute the residual for a symmetric matrix times
            // its inverse.

            dpot03(UPLO, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));
            NT = 2;
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

          // Do only the condition estimate if INFO.value is not 0.

          if (TRFCON) {
            RCONDC.value = ZERO;
            break;
          }

          // Do for each value of NRHS in NSVAL.

          for (var IRHS = 1; IRHS <= NNS; IRHS++) {
            final NRHS = NSVAL[IRHS];

            // +    TEST 3 ( Using TRS)
            // Solve and compute residual for  A * X = B.

            // Choose a set of NRHS random solution vectors
            // stored in XACT and set up the right hand side B

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'DSYTRS';
            dsytrs(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(),
                LDA, INFO);

            // Check error code from DSYTRS and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSYTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

            // Compute the residual for the solution

            dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RESULT(3));

            // +    TEST 4 (Using TRS2)

            // Solve and compute residual for  A * X = B.

            // Choose a set of NRHS random solution vectors
            // stored in XACT and set up the right hand side B

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'DSYTRS2';
            dsytrs2(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(),
                LDA, WORK, INFO);

            // Check error code from DSYTRS2 and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSYTRS2', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

            // Compute the residual for the solution

            dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RESULT(4));

            // +    TEST 5
            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(5));

            // +    TESTS 6, 7, and 8
            // Use iterative refinement to improve the solution.

            srnamc.SRNAMT = 'DSYRFS';
            dsyrfs(
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
                IWORK(N + 1),
                INFO);

            // Check error code from DSYRFS and handle error.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSYRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(6));
            dpot05(
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

          // +    TEST 9
          // Get an estimate of RCOND.value = 1/CNDNUM.

          // }
          final ANORM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
          srnamc.SRNAMT = 'DSYCON';
          dsycon(UPLO, N, AFAC.asMatrix(), LDA, IWORK, ANORM, RCOND, WORK,
              IWORK(N + 1), INFO);

          // Check error code from DSYCON and handle error.

          if (INFO.value != 0) {
            alaerh(PATH, 'DSYCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          // Compute the test ratio to compare values of RCOND.value

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
