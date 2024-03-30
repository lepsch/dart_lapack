import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dspcon.dart';
import 'package:lapack/src/dsprfs.dart';
import 'package:lapack/src/dsptrf.dart';
import 'package:lapack/src/dsptri.dart';
import 'package:lapack/src/dsptrs.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
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
import 'dppt02.dart';
import 'dppt03.dart';
import 'dppt05.dart';
import 'dspt01.dart';

void dchksp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
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
  const NTYPES = 10;
  const NTESTS = 8;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}SP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrsy(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];
        final PACKIT = lsame(UPLO, 'U') ? 'C' : 'R';

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
            dlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            PACKIT, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from DLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // For types 3-6, zero one or more rows and columns of
        // the matrix to test that INFO is returned correctly.

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

        // Compute the L*D*L' or U*D*U' factorization of the matrix.

        final NPP = N * (N + 1) ~/ 2;
        dcopy(NPP, A, 1, AFAC, 1);
        srnamc.SRNAMT = 'DSPTRF';
        dsptrf(UPLO, N, AFAC, IWORK, INFO);

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

        // Check error code from DSPTRF.

        if (INFO.value != K) {
          alaerh(PATH, 'DSPTRF', INFO.value, K, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }
        final TRFCON = INFO.value != 0;

        // +    TEST 1
        // Reconstruct matrix from factors and compute residual.

        dspt01(UPLO, N, A, AFAC, IWORK, AINV.asMatrix(), LDA, RWORK, RESULT(1));

        // +    TEST 2
        // Form the inverse and compute the residual.

        final int NT;
        final RCONDC = Box(0.0);
        if (!TRFCON) {
          dcopy(NPP, AFAC, 1, AINV, 1);
          srnamc.SRNAMT = 'DSPTRI';
          dsptri(UPLO, N, AINV, IWORK, WORK, INFO);

          // Check error code from DSPTRI.

          if (INFO.value != 0) {
            alaerh(PATH, 'DSPTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          dppt03(
              UPLO, N, A, AINV, WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));
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
                ' UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
            NFAIL++;
          }
        }
        NRUN += NT;

        // Do only the condition estimate if INFO is not 0.

        if (TRFCON) {
          RCONDC.value = ZERO;
        } else {
          for (var IRHS = 1; IRHS <= NNS; IRHS++) {
            final NRHS = NSVAL[IRHS];

            // +    TEST 3
            // Solve and compute residual for  A * X = B.

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'DSPTRS';
            dsptrs(UPLO, N, NRHS, AFAC, IWORK, X.asMatrix(), LDA, INFO);

            // Check error code from DSPTRS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSPTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                RWORK, RESULT(3));

            // +    TEST 4
            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(4));

            // +    TESTS 5, 6, and 7
            // Use iterative refinement to improve the solution.

            srnamc.SRNAMT = 'DSPRFS';
            dsprfs(
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
                RWORK,
                RWORK(NRHS + 1),
                WORK,
                IWORK(N + 1),
                INFO);

            // Check error code from DSPRFS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DSPRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(5));
            dppt05(UPLO, N, NRHS, A, B.asMatrix(), LDA, X.asMatrix(), LDA,
                XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(6));

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = 3; K <= 7; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += 5;
          }
        }

        // +    TEST 8
        // Get an estimate of RCOND = 1/CNDNUM.

        final ANORM2 = dlansp('1', UPLO, N, A, RWORK);
        final RCOND = Box(0.0);
        srnamc.SRNAMT = 'DSPCON';
        dspcon(UPLO, N, AFAC, IWORK, ANORM2, RCOND, WORK, IWORK(N + 1), INFO);

        // Check error code from DSPCON.

        if (INFO.value != 0) {
          alaerh(PATH, 'DSPCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }

        RESULT[8] = dget06(RCOND.value, RCONDC.value);

        // Print the test ratio if it is >= THRESH.

        if (RESULT[8] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.println(
              ' UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${8.i2}, ratio =${RESULT[8].g12_5}');
          NFAIL++;
        }
        NRUN++;
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
